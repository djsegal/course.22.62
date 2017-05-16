
clear
clc
% INPUT PARAMETERS
R = 5.7; % major radius
r = 1.6; % minor radius
k = 2; % plasma elongation
P_fusion = 500; % comes from fusion power in [MW]
Cold_mass = 10135; % comes from total volume of magnets x density (~8000 kg/m3) in [tons] - CHANGES EACH RUN
N_cl = 26; % number of current leads (constant)


% HEAT LOADS
Q_15 = (14.4*(P_fusion/500)) + (35.18*(Cold_mass/10135)) + (12.4*(Cold_mass/10135)) + ... 
    (6.4*(R*r^2/24.8)) + (11.4*((14.4*(P_fusion/500))+(35.18*(Cold_mass/10135)))/(14.4+35.18)) +(4*(N_cl/30)); % total hydrogen cooling @ 15K [kW]
Q_eq1 = Q_15/3.47; % equivalen cooling power at 4.5K [kW]

Q_80 = 800*((1+k^2)*R*r)/48.2; % heating due to thermal shielding @ 80K [kW]
Q_eq2 = Q_80/24; % equivalent cooling power at 4.5K [kW]

Q_total = Q_eq1 + Q_eq2; % total cooling equivalent at 4.5 K to use in Green formula [kW]

% COSTING
Cryostat_Cost = (75.8*10^6)*((1+k^2)*R*r)/48.2; % cost of cryostat scaled to ITER
Cryoplant_Cost = (2.6*10^6*(Q_total)^0.63); % cost of cryoplant scaled to various heat loads, scaled to 4.5K and estimated using Green's Formula

Cryo_Operation = (Q_total/(15/(300-15)))/(0.141*Q_total^0.26)*(365*24)*(0.1)*1.08; % cost of power @ 0.1 $/kWhr to power the cryoplant annually
