clear all 
close all
clc


lifetime = 50;  %in years

costsbyyear = zeros(4,lifetime); %rows for: CAPEX, int, O&M, Fuel
tcostbyyear = sum(costsbyyear,1);

kWhrbyyear = zeros(1,lifetime);

r = .08;    %discount rate
discountfactor = zeros(1,lifetime);
for i = 1:lifetime
    discountfactor(i) = (1+r)^i;        %this is a matrix 1xlifetime in size that tells you (1+r)^t each year
end

LCoE = sum(tcostsbyyear./discountfactor)/sum(kWhrbyyear./discountfactor));
