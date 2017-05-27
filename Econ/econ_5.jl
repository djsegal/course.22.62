using DataArrays, DataFrames, Grid, Plots
include("divertor_material_cost.jl")
include("blah_CS")

# stuff to comment out once this code can see all the other variables in the master code
P_th = 1  # in MW
P_F = 1
R_0 = 3.3
a() = 1
kappa = 1.8
blanket_thickness() = .7
# END stuff to comment out once this code can see all the other variables in the master code

eta_BOP = .5 # we want to be able to sweep over this, from .4 - .6
T = 50 #lifetime of plant from start of construction to end of decomissioning -- years
cT = 7 #construction time -- years
MW_e = P_th * eta_BOP  #power output of the plant in MW
discountRate = .08 #discount rate for the project
Availability = .7
efficency = .6

#%%% DEFINE THE TABLE!! %%%%%%%%%%%%%%%%%%%%%
MCT = DataFrame()
MCT[:Item] = ["Magnets", "Blanket", "First_Wall", "Be_multiplier", "RF", "BOP_turbine", "BOP_heatRej", "BOP_ElecGen", "BOP_other", "MagnetAux", "Divertor", "DivertorCool", "VacuumVessel","MachineAssem", "Cryostat","ThermalShields","VacuumPumpingandFueling","CoolingWater","TritiumPlant","Cryoplant","PowerSupplies","Buildings", "WasteTreatment", "RadiologicalProtection","Land"];
MCT[:DollarsPer] = ones(size(MCT[:Item]))
MCT[:Quantity] = ones(size(MCT[:Item]))
MCT[:FabFactor] = ones(size(MCT[:Item]))
#%%% DEFINE THE TABLE!! %%%%%%%%%%%%%%%%%%%%%
fill!(zeros(size(MCT[:Item])), 1)

#%%% SET COSTS %%%%%%%%%%%%%%%%%%%%%
MCT[MCT[:Item] .== "Magnets",2] = 1 #total cost of raw materials for magnets from Cost_Total from Magnets_2262_Final.m
MCT[MCT[:Item] .== "Magnets",4] = 2 #arbitrary default FF

MCT[MCT[:Item] .== "MagnetPower",2] = 1 #Magnet code?
MCT[MCT[:Item] .== "MagnetPower",3] = 1 #TO BE FILLED IN

MCT[MCT[:Item] .== "RF",2] = 1 # cost per kW given from RF team, don't have their code yet
MCT[MCT[:Item] .== "RF",3] = 1 # number of kW given from RF team, don't have their code yet

MCT[MCT[:Item] .== "BOP_turbine",2] = 360e6*(P_th/2000)^.8*(eta_BOP/.6)

MCT[MCT[:Item] .== "BOP_heatRej",2] =  87e6 * (1- eta_BOP) * P_th / 2300

MCT[MCT[:Item] .== "BOP_ElecGen",2] = 183e6(MW_e/1200)^.5

MCT[MCT[:Item] .== "BOP_other",2] =  197e6*(MW_e/1000)^.8

MCT[MCT[:Item] .== "Divertor",2] = divertor_material_cost()
# want this to be the output from the
#function divertor_material_cost which is in the divertor folder of the Final Code Blocks folder
#the above function replaces cost.m that Will sent
MCT[MCT[:Item] .== "Divertor",4] = 5 #arbitrary default FF

MCT[MCT[:Item] .== "DivertorCool",2] = 50e6 #cst for any divertor

MCT[MCT[:Item] .== "VacuumPumpingandFueling",2] = 78e6*(P_F/1758)^.8+78e6 #from aries

MCT[MCT[:Item] .== "MachineAssem",2] = 82.6e3 * 1.34 * 1552 * R_0 / 6.2 * 1.10 #1.10 is the inflation conversion

MCT[MCT[:Item] .== "Cryostat",2] = blah_CS() #"Cryostat_Cost()"  #from magnet code
MCT[MCT[:Item] .== "Cryostat",4] = 1.2 # arbitrary default FF

#assume that FW, Be, VV, Blanket, and Thermal Shield are hollow elliptic tori
b = a()*kappa
MCT[MCT[:Item] .== "First_Wall",2] = 29 * 19250 # $$/kg * kg/m^3 of tungsten
MCT[MCT[:Item] .== "First_Wall",3] = 2*pi^2* R_0* ( (a()+.005)*(b+.005) - a()*b ) ##volume
MCT[MCT[:Item] .== "First_Wall",4] = 10 #arbitrary default FF

MCT[MCT[:Item] .== "Be_multiplier",2] = 257 * 1850  # $$/kg * kg/m^3 of beryllium
MCT[MCT[:Item] .== "Be_multiplier",3] = 2*pi^2* R_0* ( (a()+.0804)*(b+0.0804) - (a()+0.0304)*(b+0.0304) ) #volume
MCT[MCT[:Item] .== "Be_multiplier",4] = 10 #arbitrary default FF

MCT[MCT[:Item] .== "VacuumVessel",2] = 56 * 8192  #T$$/kg * kg/m^3 of inconel 718
MCT[MCT[:Item] .== "VacuumVessel",3] = 2*pi^2* R_0* ( (a()+.0304)*(b+.0304) - (a()+.005)*(b+.005) + (a()+0.7258)*(b+0.7258) - (a()+0.7004)*(b+0.7004) ) #volume
MCT[MCT[:Item] .== "VacuumVessel",4] = 20 #arbitrary default FF

#multiply volume of flibe needed for tank by two to account for HX etc
MCT[MCT[:Item] .== "Blanket",2] = 154 * 1940 # T$$/kg * kg/m^3 of Flibe
MCT[MCT[:Item] .== "Blanket",3] = 2* 2*pi^2* R_0* ( (a()+0.7004)*(b+0.7004) - (a()+0.0804)*(b+0.0804) ) #volume
MCT[MCT[:Item] .== "Blanket",4] = 1.8 #min set by arc for tank etc

MCT[MCT[:Item] .== "ThermalShields",2] = 26.4 * 3760 # T$$/kg * kg/m^3 of TiH2
MCT[MCT[:Item] .== "ThermalShields",3] = 2*pi^2* R_0* ( (a()+0.7258 + shield_thickness())*(b+0.7258 + shield_thickness()) - (a()+0.7258)*(b+0.7258) ) #volume
MCT[MCT[:Item] .== "ThermalShields",4] = 5 #arbitrary default FF

MCT[MCT[:Item] .== "TritiumPlant",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "TritiumPlant",3] = 1 #TO BE FILLED IN

MCT[MCT[:Item] .== "Cryoplant",2] = blah_CS() #"Cryoplant_Cost()"  #from magnet code
MCT[MCT[:Item] .== "Cryoplant",4] = 1.2 #arbitrary default FF

# MCT[MCT[:Item] .== "PowerSupplies",2] = 4200e3
# MCT[MCT[:Item] .== "PowerSupplies",3] = MW_e
# MCT[MCT[:Item] .== "PowerSupplies",4] = 1

MCT[MCT[:Item] .== "Buildings",2] = 1.292e6 #
MCT[MCT[:Item] .== "Buildings",3] = MW_e #

MCT[MCT[:Item] .== "RadiologicalProtection",2] = 1.85e6

MCT[MCT[:Item] .== "Land",2] = 90e3 #
MCT[MCT[:Item] .== "Land",3] = MW_e #

#%%% SET COSTS %%%%%%%%%%%%%%%%%%%%%

#%%% CALCULATE THE Total Cost per Item %%%%%%%%%%%%%%%%%%%%%
MCT[:TotalCost] = MCT[:DollarsPer] .* MCT[:Quantity] .* MCT[:FabFactor]

#%%% CALCULATE THE COST / MW %%%%%%%%%%%%%%%%%%%%%
OvernightCost = sum(MCT[:TotalCost])
OvernightCostperMW = OvernightCost / MW_e

#%%% spread costs over construction time %%%%%%%%%%%%%%%%%%%%%
ParsonDist = [.1, .2, .4, .2, .1;]
ParsonYears = 1.0:length(ParsonDist)
dummyYears = linspace(1,length(ParsonDist),cT)
ig = InterpGrid( ParsonDist, BCnil, InterpLinear)
y=[ig[x] for x in dummyYears]
ConstructDist = y/sum(y)

#%%% define cost schedule %%%%%%%%%%%%%%%%%%%%%
CS = zeros(3, T + cT)
# row 1 is capital cost, row 2 is operations, and row 3 is decomissioning

#%%% add in construction costs %%%%%%%%%%%%%%%%%%%%%
CS[1,1:cT] = ConstructDist * OvernightCost
kW_h_per_year = Availability*MW_e*1000*3.15e7/3600
CS[2,cT+1:end] = 540e3*MW_e + 3.5e-4*kW_h_per_year  +  blah_CS() #blah_CS should be Cryo_Operation()
#decommissioning cost is scaled with lifetime and converted to 2017 dollars from 2009 $
# WE ARE NEGLECTING FUEL COSTS! THATS FINE THOUGH -- 600*P_F ?

#%%% define total costs PER YEAR!!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
TCy = zeros(2, T + cT)
TCy[1,:] = sum(CS,1)

#%%% define electricity produced per year in MWh %%%%%%%%%%%%%%%%%%%%%
ElecProd = zeros(2, T + cT)
ElecProd[1,cT+1:end] = fill!(zeros(1, T), Availability * MW_e * 365 * 24)


# %%% discount cost per year and electricity per year %%%%%%%%%%%%%%%%%%%%%
for i = 1:(cT + T)
  ElecProd[2,i] = ElecProd[1,i]/(1.+discountRate)^i
  TCy[2,i] = TCy[1,i]/(1.+discountRate)^i
end

########## CALCULATE LEVELIZED COST OF ELECTRICITY #########
ElecPV = sum(ElecProd[2,:])
CostPV = sum(TCy[2,:]);
lcoe = CostPV / ElecPV; # units are dollars per MW-hr
print(lcoe)
