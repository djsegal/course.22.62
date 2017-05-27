Pkg.add("PyPlot")
using PyPlot

using Polynomials

using QuadGK

using NLsolve

using Elliptic

using Unitful
using Unitful.DefaultSymbols

using Tokamak

Tokamak.load_input(" R_0 = 3.3 * 1u\"m\" ")
Tokamak.load_input(" B_0 = 9.2 * 1u\"T\" ")
Tokamak.load_input(" I_M = 8 * 1u\"MA\" ")
Tokamak.load_input(" T_k = 15 * 1u\"keV\" ")
Tokamak.load_input(" n_bar = 1.5 * 1u\"n20\" ")
Tokamak.load_input(" delta = 0.45 ")
Tokamak.load_input(" epsilon = 0.3424242424 ")
Tokamak.load_input(" enable_blanket_derive = false ")


###############################################
###############################################
## CENTRAL SOLENOID Dimensions and Requirements
## Parameters

## Solenoid Thicknesses and Volumes
cMCS = da - hts_thickness() # Thickness of Structure
a2CS = a1+da # Outer Radius
VJ_CS = pi*( ((a2CS+a1)/2 + hts_thickness()/2)^2 - ((a2CS+a1)/2 - hts_thickness()/2)^2)*Tokamak.solenoid_length() # Volume of HTS
VM_CS = pi*( (a2CS)^2 - ((a2CS+a1)/2 + hts_thickness()/2)^2 + ((a2CS+a1)/2 - hts_thickness()/2)^2 - a1^2)*Tokamak.solenoid_length() # Volume of Structure


Vol_CS = VJ_CS + VM_CS
Cost_CS = Tokamak.Price_St*(VM_CS+(1-Tokamak.magnet_hts_fraction)*VJ_CS)*8000 + Tokamak.magnet_hts_fraction*VJ_CS*(Tokamak.Price_HTS/Tokamak.Area_Tape())

Vol_ST_Total = (VM_CS + (sum(Vst_PF())) + (sum(Vsc_PF())*(1-Tokamak.magnet_hts_fraction)) + VJ_CS*(1-Tokamak.magnet_hts_fraction) + Tokamak.V_TF() + Tokamak.Vol_WP()*(1-Tokamak.Frac_HTS()))
Cost_ST_Total = Vol_ST_Total*Tokamak.Price_St*8000

Vol_HTS_Total = (VJ_CS*Tokamak.magnet_hts_fraction+(sum(Vsc_PF())*Tokamak.magnet_hts_fraction)+Tokamak.Vol_WP()*Tokamak.Frac_HTS())
Cost_HTS_Total = Vol_HTS_Total*Tokamak.Price_HTS/(Tokamak.Area_Tape())

Cost_Total = Cost_HTS_Total+Cost_ST_Total
Cost_Total11 = Cost_HTS_Total + Cost_ST_Total
