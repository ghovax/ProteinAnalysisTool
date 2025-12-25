//! PyMOL-style selection language parser

use super::SelectionExpression;

#[derive(Debug, Clone, PartialEq)]
enum Token {
    ChainKeyword,
    ResidueNameKeyword,
    ResidueIdentifierKeyword,
    AtomNameKeyword,
    ElementKeyword,
    BackboneKeyword,
    SidechainKeyword,
    HelixKeyword,
    SheetKeyword,
    WithinKeyword,
    OfKeyword,
    AllKeyword,
    AndOperator,
    OrOperator,
    NotOperator,
    LeftParenthesis,
    RightParenthesis,
    Identifier(String),
    Number(f32),
    ResidueRange(isize, isize),
}

pub struct SelectionParser<'a> {
    token_list: Vec<Token>,
    current_position: usize,
    _input_string: &'a str,
}

impl<'a> SelectionParser<'a> {
    pub fn new(input_string: &'a str) -> Self {
        let token_list = tokenize_selection_string(input_string);
        Self {
            token_list,
            current_position: 0,
            _input_string: input_string,
        }
    }

    pub fn parse_expression(&mut self) -> Result<SelectionExpression, String> {
        if self.token_list.is_empty() {
            return Ok(SelectionExpression::All);
        }
        let expression = self.parse_or_expression()?;
        if self.current_position < self.token_list.len() {
            return Err(format!("Unexpected token at end of expression: {:?}", self.token_list[self.current_position]));
        }
        Ok(expression)
    }

    fn parse_or_expression(&mut self) -> Result<SelectionExpression, String> {
        let mut left_expression = self.parse_and_expression()?;
        while self.peek_next_token() == Some(&Token::OrOperator) {
            self.advance_token_position();
            let right_expression = self.parse_and_expression()?;
            left_expression = SelectionExpression::Or(Box::new(left_expression), Box::new(right_expression));
        }
        Ok(left_expression)
    }

    fn parse_and_expression(&mut self) -> Result<SelectionExpression, String> {
        let mut left_expression = self.parse_not_expression()?;
        // Support both explicit 'and' and implicit 'and' (space-separated)
        while let Some(next_token) = self.peek_next_token() {
            if next_token == &Token::AndOperator {
                self.advance_token_position();
                let right_expression = self.parse_not_expression()?;
                left_expression = SelectionExpression::And(Box::new(left_expression), Box::new(right_expression));
            } else if matches!(next_token, Token::ChainKeyword | Token::ResidueNameKeyword | Token::ResidueIdentifierKeyword | Token::AtomNameKeyword | Token::ElementKeyword | 
                               Token::BackboneKeyword | Token::SidechainKeyword | Token::HelixKeyword | Token::SheetKeyword | 
                               Token::WithinKeyword | Token::AllKeyword | Token::NotOperator | Token::LeftParenthesis | Token::Identifier(_)) {
                let right_expression = self.parse_not_expression()?;
                left_expression = SelectionExpression::And(Box::new(left_expression), Box::new(right_expression));
            } else {
                break;
            }
        }
        Ok(left_expression)
    }

    fn parse_not_expression(&mut self) -> Result<SelectionExpression, String> {
        if self.peek_next_token() == Some(&Token::NotOperator) {
            self.advance_token_position();
            let inner_expression = self.parse_atomic_expression()?;
            Ok(SelectionExpression::Not(Box::new(inner_expression)))
        } else {
            self.parse_atomic_expression()
        }
    }

    fn parse_atomic_expression(&mut self) -> Result<SelectionExpression, String> {
        let next_token = self.peek_next_token().cloned();
        match next_token {
            Some(Token::LeftParenthesis) => {
                self.advance_token_position();
                let inner_expression = self.parse_or_expression()?;
                if self.peek_next_token() != Some(&Token::RightParenthesis) {
                    return Err("Expected matching ')'".to_string());
                }
                self.advance_token_position();
                Ok(inner_expression)
            }
            Some(Token::ChainKeyword) => {
                self.advance_token_position();
                if let Some(Token::Identifier(chain_id)) = self.advance_token_position() {
                    Ok(SelectionExpression::Chain(chain_id))
                } else {
                    Err("Expected chain identifier after 'chain' keyword".to_string())
                }
            }
            Some(Token::ResidueNameKeyword) => {
                self.advance_token_position();
                if let Some(Token::Identifier(residue_name)) = self.advance_token_position() {
                    Ok(SelectionExpression::ResidueName(residue_name.to_uppercase()))
                } else {
                    Err("Expected residue name after 'resn' keyword".to_string())
                }
            }
            Some(Token::ResidueIdentifierKeyword) => {
                self.advance_token_position();
                match self.advance_token_position() {
                    Some(Token::Number(number)) => Ok(SelectionExpression::ResidueNumber(number as isize)),
                    Some(Token::ResidueRange(start_number, end_number)) => Ok(SelectionExpression::ResidueRange(start_number, end_number)),
                    _ => Err("Expected residue number or range after 'resi' keyword".to_string()),
                }
            }
            Some(Token::AtomNameKeyword) => {
                self.advance_token_position();
                if let Some(Token::Identifier(atom_name)) = self.advance_token_position() {
                    Ok(SelectionExpression::AtomName(atom_name.to_uppercase()))
                } else {
                    Err("Expected atom name after 'name' keyword".to_string())
                }
            }
            Some(Token::ElementKeyword) => {
                self.advance_token_position();
                if let Some(Token::Identifier(element_symbol)) = self.advance_token_position() {
                    Ok(SelectionExpression::Element(element_symbol.to_uppercase()))
                } else {
                    Err("Expected element symbol after 'elem' keyword".to_string())
                }
            }
            Some(Token::BackboneKeyword) => {
                self.advance_token_position();
                Ok(SelectionExpression::Backbone)
            }
            Some(Token::SidechainKeyword) => {
                self.advance_token_position();
                Ok(SelectionExpression::Sidechain)
            }
            Some(Token::HelixKeyword) => {
                self.advance_token_position();
                Ok(SelectionExpression::Helix)
            }
            Some(Token::SheetKeyword) => {
                self.advance_token_position();
                Ok(SelectionExpression::Sheet)
            }
            Some(Token::AllKeyword) => {
                self.advance_token_position();
                Ok(SelectionExpression::All)
            }
            Some(Token::WithinKeyword) => {
                self.advance_token_position();
                let distance_threshold = if let Some(Token::Number(number)) = self.advance_token_position() {
                    number
                } else {
                    return Err("Expected numerical distance after 'within' keyword".to_string());
                };
                if self.peek_next_token() != Some(&Token::OfKeyword) {
                    return Err("Expected 'of' keyword after distance in 'within' expression".to_string());
                }
                self.advance_token_position();
                let target_expression = self.parse_atomic_expression()?;
                Ok(SelectionExpression::Within(distance_threshold, Box::new(target_expression)))
            }
            Some(Token::Identifier(identifier_string)) => {
                self.advance_token_position();
                match identifier_string.to_lowercase().as_str() {
                    "all" => Ok(SelectionExpression::All),
                    "none" => Ok(SelectionExpression::None),
                    "backbone" => Ok(SelectionExpression::Backbone),
                    "sidechain" => Ok(SelectionExpression::Sidechain),
                    "helix" => Ok(SelectionExpression::Helix),
                    "sheet" => Ok(SelectionExpression::Sheet),
                    _ => Ok(SelectionExpression::ResidueName(identifier_string.to_uppercase())),
                }
            }
            _ => Err(format!("Unexpected token encountered: {:?}", self.peek_next_token())),
        }
    }

    fn peek_next_token(&self) -> Option<&Token> {
        self.token_list.get(self.current_position)
    }

    fn advance_token_position(&mut self) -> Option<Token> {
        if self.current_position < self.token_list.len() {
            let token = self.token_list[self.current_position].clone();
            self.current_position += 1;
            Some(token)
        } else {
            None
        }
    }
}

fn tokenize_selection_string(input_string: &str) -> Vec<Token> {
    let mut tokens = Vec::new();
    let mut characters = input_string.chars().peekable();

    while let Some(&current_char) = characters.peek() {
        match current_char {
            '(' => {
                tokens.push(Token::LeftParenthesis);
                characters.next();
            }
            ')' => {
                tokens.push(Token::RightParenthesis);
                characters.next();
            }
            '-' if is_next_character_digit(&mut characters) => {
                let numeric_string = read_while_predicate_true(&mut characters, |c| c.is_digit(10) || c == '-' || c == '.');
                if let Some(dash_index) = numeric_string.find('-').filter(|&index| index > 0) {
                    let start_number = numeric_string[..dash_index].parse().unwrap_or(0);
                    let end_number = numeric_string[dash_index+1..].parse().unwrap_or(0);
                    tokens.push(Token::ResidueRange(start_number, end_number));
                } else if let Ok(parsed_number) = numeric_string.parse::<f32>() {
                    tokens.push(Token::Number(parsed_number));
                }
            }
            c if c.is_digit(10) => {
                let numeric_string = read_while_predicate_true(&mut characters, |c| c.is_digit(10) || c == '-' || c == '.');
                if let Some(dash_index) = numeric_string.find('-') {
                    let start_number = numeric_string[..dash_index].parse().unwrap_or(0);
                    let end_number = numeric_string[dash_index+1..].parse().unwrap_or(0);
                    tokens.push(Token::ResidueRange(start_number, end_number));
                } else if let Ok(parsed_number) = numeric_string.parse::<f32>() {
                    tokens.push(Token::Number(parsed_number));
                }
            }
            c if c.is_alphabetic() => {
                let alphanumeric_string = read_while_predicate_true(&mut characters, |c| c.is_alphanumeric() || c == '_');
                match alphanumeric_string.to_lowercase().as_str() {
                    "chain" => tokens.push(Token::ChainKeyword),
                    "resn" => tokens.push(Token::ResidueNameKeyword),
                    "resi" => tokens.push(Token::ResidueIdentifierKeyword),
                    "name" => tokens.push(Token::AtomNameKeyword),
                    "elem" => tokens.push(Token::ElementKeyword),
                    "backbone" => tokens.push(Token::BackboneKeyword),
                    "sidechain" => tokens.push(Token::SidechainKeyword),
                    "helix" => tokens.push(Token::HelixKeyword),
                    "sheet" => tokens.push(Token::SheetKeyword),
                    "within" => tokens.push(Token::WithinKeyword),
                    "of" => tokens.push(Token::OfKeyword),
                    "all" => tokens.push(Token::AllKeyword),
                    "and" => tokens.push(Token::AndOperator),
                    "or" => tokens.push(Token::OrOperator),
                    "not" => tokens.push(Token::NotOperator),
                    _ => tokens.push(Token::Identifier(alphanumeric_string)),
                }
            }
            _ if current_char.is_whitespace() => {
                characters.next();
            }
            _ => {
                characters.next();
            }
        }
    }
    tokens
}

fn is_next_character_digit(characters: &mut std::iter::Peekable<std::str::Chars>) -> bool {
    let mut cloned_characters = characters.clone();
    cloned_characters.next(); // skip current
    cloned_characters.peek().map_or(false, |c| c.is_digit(10))
}

fn read_while_predicate_true<F>(characters: &mut std::iter::Peekable<std::str::Chars>, mut predicate: F) -> String 
where F: FnMut(char) -> bool {
    let mut accumulated_string = String::new();
    while let Some(&current_char) = characters.peek() {
        if predicate(current_char) {
            accumulated_string.push(current_char);
            characters.next();
        } else {
            break;
        }
    }
    accumulated_string
}

pub fn parse_selection_expression(input_string: &str) -> Result<SelectionExpression, String> {
    let mut selection_parser = SelectionParser::new(input_string);
    selection_parser.parse_expression()
}