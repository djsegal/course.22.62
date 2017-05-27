using DataArrays, DataFrames, Grid

# MW = 500;   #power output of the plant from some other group -- MW
# P_f = 1;

#%%% DEFINE THE TABLE!! %%%%%%%%%%%%%%%%%%%%%
MCT = DataFrame();
MCT[:Item] = ["Magnets", "Blanket", "RF", "BOP_turbine", "BOP_heatRej", "BOP_ElecGen", "BOP_other", "MagnetAux", "Divertor", "DivertorAux", "VacuumVessel","MachineAssem", "Cryostat","ThermalShields","VacuumPumpingandFueling","CoolingWater","TritiumPlant","Cryoplant","PowerSupplies","Buildings", "WasteTreatment", "RadiologicalProtection","Land", "ICHCD","ECHCD","NBHCD"];
MCT[:DollarsPer] = ones(size(MCT[:Item]))
MCT[:Quantity] = ones(size(MCT[:Item]));
MCT[:FabFactor] = ones(size(MCT[:Item]))
#%%% DEFINE THE TABLE!! %%%%%%%%%%%%%%%%%%%%%
fill!(zeros(size(MCT[:Item])), 1)

#%%% SET COSTS %%%%%%%%%%%%%%%%%%%%%
MCT[MCT[:Item] .== "Magnets",2] = 1 #see magnet code
MCT[MCT[:Item] .== "Magnets",3] = 1 #see magnet code
MCT[MCT[:Item] .== "Magnets",4] = 10 #arbitrary for now
MCT[MCT[:Item] .== "MagnetAux",2] = 1 #Magnet code?
MCT[MCT[:Item] .== "MagnetAux",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "MagnetAux",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "Blanket",2] = 1 #needs to be filled in
MCT[MCT[:Item] .== "Blanket",3] = 1 #needs to be filled in
MCT[MCT[:Item] .== "Blanket",4] = 100 #???
MCT[MCT[:Item] .== "RF",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "RF",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "RF",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "BOP_turbine",2] =  360000000*(Pth gross/2000)^.8*(efficency/.6) #as is
MCT[MCT[:Item] .== "BOP_turbine",3] = 1 #number of turbines i guess?
MCT[MCT[:Item] .== "BOP_turbine",4] = 1 #Placeholder
MCT[MCT[:Item] .== "BOP_heatRej",2] =  67000000 #as is
MCT[MCT[:Item] .== "BOP_heatRej",3] = 1 #number of turbines i guess?
MCT[MCT[:Item] .== "BOP_heatRej",4] = 1 #Placeholder
MCT[MCT[:Item] .== "BOP_ElecGen",2] = 183000000(Pe/1200)^.5 #as is
MCT[MCT[:Item] .== "BOP_ElecGen",3] = 1
MCT[MCT[:Item] .== "BOP_ElecGen",4] = 1 #Placeholder
MCT[MCT[:Item] .== "BOP_other",2] =  197000000*(Pe/1000)^.8 #as is
MCT[MCT[:Item] .== "BOP_other",3] = 1
MCT[MCT[:Item] .== "BOP_other",4] = 1 #Placeholder
MCT[MCT[:Item] .== "Divertor",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "Divertor",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "Divertor",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "DivertorAux",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "DivertorAux",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "DivertorAux",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "VacuumPumpingandFueling",2] = 78000000*(P_f/1758)^.8+7000000
MCT[MCT[:Item] .== "VacuumPumpingandFueling",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "VacuumPumpingandFueling",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "MachineAssem",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "MachineAssem",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "MachineAssem",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "Cryostat",2] = 1 #from magnet code
MCT[MCT[:Item] .== "Cryostat",3] = 1 #from magnet code
MCT[MCT[:Item] .== "Cryostat",4] = 1 #from magnet code
MCT[MCT[:Item] .== "ThermalShields",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "ThermalShields",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "ThermalShields",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "VacuumVessel",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "VacuumVessel",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "VacuumVessel",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "CoolingWater",2] = (125000000*P_th/2000)^0.55
MCT[MCT[:Item] .== "CoolingWater",3] = 1
MCT[MCT[:Item] .== "CoolingWater",4] = 1 #arbitrary
MCT[MCT[:Item] .== "TritiumPlant",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "TritiumPlant",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "TritiumPlant",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "Cryoplant",2] = 1 #from magnet code
MCT[MCT[:Item] .== "Cryoplant",3] = 1 #from magnet code
MCT[MCT[:Item] .== "Cryoplant",4] = 1 #from magnet code
MCT[MCT[:Item] .== "PowerSupplies",2] = 4200000
MCT[MCT[:Item] .== "PowerSupplies",3] = MW
MCT[MCT[:Item] .== "PowerSupplies",4] = 1
MCT[MCT[:Item] .== "Buildings",2] = 1292000 #
MCT[MCT[:Item] .== "Buildings",3] = MW #
MCT[MCT[:Item] .== "Buildings",4] = 1 #
MCT[MCT[:Item] .== "WasteTreatment",2] = 0
MCT[MCT[:Item] .== "WasteTreatment",3] = 1
MCT[MCT[:Item] .== "WasteTreatment",4] = 1
MCT[MCT[:Item] .== "RadiologicalProtection",2] = 1850000
MCT[MCT[:Item] .== "RadiologicalProtection",3] = 1
MCT[MCT[:Item] .== "RadiologicalProtection",2] = 1
MCT[MCT[:Item] .== "Land",2] = 90000 #
MCT[MCT[:Item] .== "Land",3] = MW #
MCT[MCT[:Item] .== "Land",4] = 1 #
MCT[MCT[:Item] .== "ICHCD",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "ICHCD",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "ICHCD",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "ECHCD",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "ECHCD",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "ECHCD",4] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "NBHCD",2] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "NBHCD",3] = 1 #TO BE FILLED IN
MCT[MCT[:Item] .== "NBHCD",4] = 1 #TO BE FILLED IN
#%%% SET COSTS %%%%%%%%%%%%%%%%%%%%%


#%%% CALCULATE THE Total Cost per Item %%%%%%%%%%%%%%%%%%%%%
MCT[:TotalCost] = MCT[:DollarsPer] .* MCT[:Quantity] .* MCT[:FabFactor]

#%%% CALCULATE THE COST / MW %%%%%%%%%%%%%%%%%%%%%
OvernightCost = sum(MCT[:TotalCost])
OvernightCostperMW = OvernightCost / MW;


#%%% spread costs over construction time %%%%%%%%%%%%%%%%%%%%%
ParsonDist = [.1, .2, .4, .2, .1;]
ParsonYears = 1.0:length(ParsonDist);
dummyYears = linspace(1,length(ParsonDist),construction_time_in_years);
ig = InterpGrid( ParsonDist, BCnil, InterpLinear);
y=[ig[x] for x in dummyYears];
ConstructDist = y/sum(y);

#%%% define cost schedule %%%%%%%%%%%%%%%%%%%%%
CS = zeros(3, plant_life_in_years + construction_time_in_years)

#%%% add in construction costs %%%%%%%%%%%%%%%%%%%%%
CS[1,1:construction_time_in_years] = ConstructDist * OvernightCost;
CS[2,construction_time_in_years+1:end] = 540000*MW + 600*P_f


#%%% define total costs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
TC = sum(CS,1)
# CS[4,:] = TC;
# TC[1,1:plant_life_in_years] = ConstructionDist * OvernightCost + 540000*MW + 600*Pf

#%%% define electricity produced per year in MWh %%%%%%%%%%%%%%%%%%%%%
ElecProd = zeros(2, plant_life_in_years + construction_time_in_years)
ElecProd[1,construction_time_in_years+1:end] = fill!(zeros(1, plant_life_in_years), econ_availability * MW * 365 * 24)


#%%% discount cost per year and electricity per year %%%%%%%%%%%%%%%%%%%%%
for i = 1:construction_time_in_years + plant_life_in_years
  ElecProd[2,i] = ElecProd[1,i]/(1.+discount_rate)^i
  TC[2,i] = sum(CS[1:3,i])/(1.+discount_rate)^i
end


TotalYears = linspace(1,construction_time_in_years + plant_life_in_years, construction_time_in_years + plant_life_in_years);
ElecPV = sum(ElecProd[2,:])
CostPV = sum(TC[2,:]);
lcoe = CostPV / ElecPV;
