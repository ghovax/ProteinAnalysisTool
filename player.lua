-- player.lua
local t = os.clock()
local speed = 3.0
local radius = 0.45

local x = math.sin(2 * t * speed) * radius
local y = math.cos(2 * t * speed) * radius

-- Print to terminal so we can see it's actually working
print("Moving to: " .. x .. ", " .. y)

transform:set_pos(x, y, 0.0)