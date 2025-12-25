//! Secondary structure annotations from HELIX and SHEET records
//!
//! This module defines structures for storing secondary structure information
//! parsed from PDB HELIX/SHEET records and mmCIF _struct_conf/_struct_sheet_range

/// The class of a helix as defined in PDB format
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HelixClass {
    /// Right-handed alpha helix (class 1, most common)
    #[default]
    RightHandedAlpha = 1,
    /// Right-handed omega helix (class 2)
    RightHandedOmega = 2,
    /// Right-handed pi helix (class 3)
    RightHandedPi = 3,
    /// Right-handed gamma helix (class 4)
    RightHandedGamma = 4,
    /// Right-handed 3-10 helix (class 5)
    RightHanded310 = 5,
    /// Left-handed alpha helix (class 6)
    LeftHandedAlpha = 6,
    /// Left-handed omega helix (class 7)
    LeftHandedOmega = 7,
    /// Left-handed gamma helix (class 8)
    LeftHandedGamma = 8,
    /// 2-7 ribbon/helix (class 9)
    Ribbon27 = 9,
    /// Polyproline (class 10)
    Polyproline = 10,
    /// Unknown helix class
    Unknown = 0,
}

impl HelixClass {
    /// Create a HelixClass from the PDB class number
    pub fn from_class_number(class_number: usize) -> Self {
        match class_number {
            1 => HelixClass::RightHandedAlpha,
            2 => HelixClass::RightHandedOmega,
            3 => HelixClass::RightHandedPi,
            4 => HelixClass::RightHandedGamma,
            5 => HelixClass::RightHanded310,
            6 => HelixClass::LeftHandedAlpha,
            7 => HelixClass::LeftHandedOmega,
            8 => HelixClass::LeftHandedGamma,
            9 => HelixClass::Ribbon27,
            10 => HelixClass::Polyproline,
            _ => HelixClass::Unknown,
        }
    }

    /// Get the PDB class number for this helix class
    pub fn class_number(&self) -> usize {
        match self {
            HelixClass::RightHandedAlpha => 1,
            HelixClass::RightHandedOmega => 2,
            HelixClass::RightHandedPi => 3,
            HelixClass::RightHandedGamma => 4,
            HelixClass::RightHanded310 => 5,
            HelixClass::LeftHandedAlpha => 6,
            HelixClass::LeftHandedOmega => 7,
            HelixClass::LeftHandedGamma => 8,
            HelixClass::Ribbon27 => 9,
            HelixClass::Polyproline => 10,
            HelixClass::Unknown => 0,
        }
    }
}

/// The sense (orientation) of a sheet strand relative to the previous strand
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SheetSense {
    /// First strand in the sheet (sense = 0)
    #[default]
    First = 0,
    /// Parallel to previous strand (sense = 1)
    Parallel = 1,
    /// Anti-parallel to previous strand (sense = -1)
    AntiParallel = -1,
}

impl SheetSense {
    /// Create a SheetSense from the PDB sense value
    pub fn from_sense_value(sense_value: isize) -> Self {
        match sense_value {
            0 => SheetSense::First,
            1 => SheetSense::Parallel,
            -1 => SheetSense::AntiParallel,
            _ => SheetSense::First,
        }
    }

    /// Get the PDB sense value for this sheet sense
    pub fn sense_value(&self) -> isize {
        match self {
            SheetSense::First => 0,
            SheetSense::Parallel => 1,
            SheetSense::AntiParallel => -1,
        }
    }
}

/// A helix secondary structure annotation from a PDB HELIX record
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, PartialEq)]
pub struct HelixAnnotation {
    /// Serial number of the helix (unique within the file)
    pub serial_number: usize,
    /// Helix identifier (e.g., "H1", "3-1")
    pub helix_id: String,
    /// Chain identifier where the helix is located
    pub chain_id: String,
    /// Starting residue name (e.g., "THR")
    pub start_residue_name: String,
    /// Starting residue sequence number
    pub start_residue_number: isize,
    /// Starting residue insertion code
    pub start_insertion_code: Option<String>,
    /// Ending residue name (e.g., "LEU")
    pub end_residue_name: String,
    /// Ending residue sequence number
    pub end_residue_number: isize,
    /// Ending residue insertion code
    pub end_insertion_code: Option<String>,
    /// Class of the helix (alpha, 3-10, etc.)
    pub helix_class: HelixClass,
    /// Length of the helix in residues
    pub length: usize,
}

impl HelixAnnotation {
    /// Create a new HelixAnnotation
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        serial_number: usize,
        helix_id: String,
        chain_id: String,
        start_residue_name: String,
        start_residue_number: isize,
        start_insertion_code: Option<String>,
        end_residue_name: String,
        end_residue_number: isize,
        end_insertion_code: Option<String>,
        helix_class: HelixClass,
        length: usize,
    ) -> Self {
        Self {
            serial_number,
            helix_id,
            chain_id,
            start_residue_name,
            start_residue_number,
            start_insertion_code,
            end_residue_name,
            end_residue_number,
            end_insertion_code,
            helix_class,
            length,
        }
    }

    /// Check if a residue is part of this helix
    pub fn contains_residue(&self, chain_id: &str, residue_number: isize) -> bool {
        self.chain_id == chain_id
            && residue_number >= self.start_residue_number
            && residue_number <= self.end_residue_number
    }
}

/// A sheet strand annotation from a PDB SHEET record
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, PartialEq)]
pub struct SheetAnnotation {
    /// Strand number within the sheet
    pub strand_number: usize,
    /// Sheet identifier (e.g., "S1")
    pub sheet_id: String,
    /// Number of strands in this sheet
    pub num_strands: usize,
    /// Chain identifier where the strand is located
    pub chain_id: String,
    /// Starting residue name (e.g., "CYS")
    pub start_residue_name: String,
    /// Starting residue sequence number
    pub start_residue_number: isize,
    /// Starting residue insertion code
    pub start_insertion_code: Option<String>,
    /// Ending residue name (e.g., "TRP")
    pub end_residue_name: String,
    /// Ending residue sequence number
    pub end_residue_number: isize,
    /// Ending residue insertion code
    pub end_insertion_code: Option<String>,
    /// Sense of this strand relative to the previous strand
    pub sense: SheetSense,
}

impl SheetAnnotation {
    /// Create a new SheetAnnotation
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        strand_number: usize,
        sheet_id: String,
        num_strands: usize,
        chain_id: String,
        start_residue_name: String,
        start_residue_number: isize,
        start_insertion_code: Option<String>,
        end_residue_name: String,
        end_residue_number: isize,
        end_insertion_code: Option<String>,
        sense: SheetSense,
    ) -> Self {
        Self {
            strand_number,
            sheet_id,
            num_strands,
            chain_id,
            start_residue_name,
            start_residue_number,
            start_insertion_code,
            end_residue_name,
            end_residue_number,
            end_insertion_code,
            sense,
        }
    }

    /// Check if a residue is part of this sheet strand
    pub fn contains_residue(&self, chain_id: &str, residue_number: isize) -> bool {
        self.chain_id == chain_id
            && residue_number >= self.start_residue_number
            && residue_number <= self.end_residue_number
    }
}
