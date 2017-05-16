%Takes a table of LCOE and plots it against a table of B values

B = 0 %code for accessing and creating the B vector
LCOE = 0 %code for accessing and creating the vector of LCOE

Fission = ones(1, length(B));
Fission = Fission*(124/1000);

Coal = ones(1, length(B));
Coal = Coal*(100/1000);

Gas_com_cycle = ones(1, length(B));
Gas_com_cycle = Gas_com_cycle*(60/1000);

Wind = ones(1, length(B));
Wind = Wind*(54/1000);

Solar = ones(1, length(B));
Solar = Solar*(43/1000);

figure(1)
hold on 
title('LCOE versus B')
xlabel('B(T)')
ylabel('LCOE(USD/kwh)')
plot(B,LCOE)
plot(B,Fission)
plot(B,Coal)
plot(B,Gas_com_cycle)
plot(B,Wind)
plot(B,Solar)
legend('Fusion','Fission','Coal','Gas Combined Cycle','Wind','Solar')
hold off



















