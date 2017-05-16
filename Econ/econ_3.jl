using DataArrays, DataFrames, Grid, Plots

T = 50; #lifetime of plant from start of construction to end of decomissioning -- years
cT = 7; #construction time -- years
MW = 500;   #power output of the plant from some other group -- MW
discountRate = .08; #discount rate for the project
Availability = .7;

#%%% DEFINE THE TABLE!! %%%%%%%%%%%%%%%%%%%%%
MCT = DataFrame();
MCT[:Item] = ["Magnets", "Blanket", "RF", "BOP_turbine", "BOP_heatRej", "BOP_ElecGen", "BOP_other", "MagnetAux", "Divertor", "DivertorAux", "VacuumVessel","MachineAssem", "Cryostat","ThermalShields","VacuumPumpingandFueling","CoolingWater","TritiumPlant","Cryoplant","PowerSupplies","Buildings", "WasteTreatment", "RadiologicalProtection","Land", "ICHCD","ECHCD","NBHCD"];
MCT[:DollarsPer] = ones(size(MCT[:Item]))
MCT[:Quantity] = ones(size(MCT[:Item]));
MCT[:FabFactor] = ones(size(MCT[:Item]))
#%%% DEFINE THE TABLE!! %%%%%%%%%%%%%%%%%%%%%
fill!(zeros(size(MCT[:Item])), 1)

#%%% SET COSTS %%%%%%%%%%%%%%%%%%%%%
MCT[MCT[:Item] .== "VacuumVessel",2] = 78000000*(P_f/1758)^.8+7000000
MCT[MCT[:Item] .== "CoolingWater",2] = (125000000*P_th/2000)^0.55
MCT[MCT[:Item] .== "WasteTreatment",2] = 0
MCT[MCT[:Item] .== "RadiologicalProtection",2] = 0
MCT[MCT[:Item] .== "Magnets",3] = 18
MCT[MCT[:Item] .== "Blanket",4] = 100
#%%% SET COSTS %%%%%%%%%%%%%%%%%%%%%


#%%% CALCULATE THE Total Cost per Item %%%%%%%%%%%%%%%%%%%%%
MCT[:TotalCost] = MCT[:DollarsPer] .* MCT[:Quantity] .* MCT[:FabFactor]

#%%% CALCULATE THE COST / MW %%%%%%%%%%%%%%%%%%%%%
OvernightCost = sum(MCT[:TotalCost])
OvernightCostperMW = OvernightCost / MW;


#%%% spread costs over construction time %%%%%%%%%%%%%%%%%%%%%
ParsonDist = [.1, .2, .4, .2, .1;]
ParsonYears = 1.0:length(ParsonDist);
dummyYears = linspace(1,length(ParsonDist),cT);
ig = InterpGrid( ParsonDist, BCnil, InterpLinear);
y=[ig[x] for x in dummyYears];
ConstructDist = y/sum(y);

#%%% define cost schedule %%%%%%%%%%%%%%%%%%%%%
CS = zeros(4, T + cT)

#%%% add in construction costs %%%%%%%%%%%%%%%%%%%%%
CS[1,1:cT] = ConstructDist * OvernightCost;


#%%% define electricity produced per year in MWh %%%%%%%%%%%%%%%%%%%%%
ElecProd = zeros(2, T + cT)
ElecProd[1,cT+1:end] = fill!(zeros(1, T), Availability * MW * 365 * 24)


#%%% discount cost per year and electricity per year %%%%%%%%%%%%%%%%%%%%%
for i = 1:cT + T
  ElecProd[2,i] = ElecProd[1,i]/(1.+discountRate)^i
  CS[4,i] = sum(CS[1:3,i])/(1.+discountRate)^i
end


TotalYears = linspace(1,cT + T, cT + T);
ElecPV = sum(ElecProd[2,:])
CostPV = sum(CS[4,:]);
lcoe = CostPV / ElecPV;
