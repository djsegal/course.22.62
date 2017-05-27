using Tokamak
using Unitful
using Unitful.DefaultSymbols

Tokamak.load_input(" R_0 = 5.7 * 1u\"m\" ")
Tokamak.load_input(" kappa = 2.0 ")
Tokamak.load_input(" epsilon = 0.280701754385 ")
Tokamak.load_input(" T_k = 15u\"keV\" ")

# chosen to get P_F = 500
Tokamak.load_input(" n_bar = 0.5522994335297159 * 1u\"n20\" ")

# INPUT PARAMETERS

# HEAT LOADS

# COSTING
Cryostat_Cost = (75.8*10^6)*((1+Tokamak.kappa^2)*( Tokamak.R_0 / 1u"m" )*( Tokamak.a() / 1u"m" ))/48.2 # cost of cryostat scaled to ITER
Cryoplant_Cost = (2.6*10^6*(Tokamak.Q_total())^0.63) # cost of cryoplant scaled to various heat load()s, scaled to 4.5K and estimated using Green's Formula

Cryo_Operation = (Tokamak.Q_total()/(15/(300-15)))/(0.141*Tokamak.Q_total()^0.26)*(365*24)*(0.1)*1.08 # cost of power @ 0.1 $/kWhr to power the cryoplant annually
