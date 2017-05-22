+"""
 +    nuc_heat_to_mag())
 +
 +Lorem ipsum dolor sit amet.
 +"""
 +function nuc_heat_to_mag()

    cur_a = 0.4318  #units of J/10cm
    cur_b = -0.1345 #units of /m
 +
    x_start = shield_thickness() /100 # shield thickness in cm

    x_end = 100;  #dist where neutrons are statistically zero

    cur_nuc_heat_to_mag = exp(x_end*cur_b) - exp(s_start*cur_b)

    cur_nuc_heat_to_mag *= cur_a * 10

    cur_nuc_heat_to_mag /= cur_b

    cur_nuc_heat_to_mag
    
 +end
