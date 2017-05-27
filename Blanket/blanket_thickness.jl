"""
    blanket_thickness
total thickness of everything between plasma and magnets, ie:
first wall, vacuum vessel, flibe blanket, be multiplier, and shielding.
"""
function blanket_thickness()

  blanket_thickness = shield_thickness() + 0.7258

end
