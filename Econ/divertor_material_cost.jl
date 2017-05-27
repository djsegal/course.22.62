"""
    divertor_material_cost
Heating due to thermal shielding @ 80K [kW].
"""
function divertor_material_cost()

epsilon = .25
R_0 = 1
  #we are assuming that minor radio is epsilon * R_0,
  #cuz we couldn't find the name of the variable for minor radius... sorry
  cur_a= epsilon * R_0

  #Volume tungsten in cooling
  V_w_c = 0.26*cur_a*R_0;
  cost_W = 600000  #USD/m^3

  #Volume Steel
  V_steel = 0.22*cur_a*R_0
  cost_steel = 16000;
  # % cost of Helium loop cooling 50 million

  #Flibe channel tungsten
  V_w_fl = 0.075*cur_a*R_0

  # %volume inconel
  V_inconel = cur_a*R_0
  cost_inconel = 66400

  div_mat_c = V_w_c*cost_W +  V_steel* cost_steel + V_w_fl*cost_W + V_inconel* cost_inconel

  div_mat_c
end
