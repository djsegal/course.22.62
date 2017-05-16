function F = rootsol(x)
%% Parameters
global R0 I Sysol Lsol B1 mu0 cJ Isol li
%R0 = 3.3;    % INPUT FIELD ON AXIS
%I = 8e6;    % INPUT PLASMA CURRENT
%Sy = 600e6;  % INPUT MAXIMUM ALLOWABLE STRESS
%Jmax = 75e6; % INPUT MAXIMUM ALLOWABLE CURRENT DENSITY
%L = 6;         % INPUT LENGTH OF SOLENOID
%B1 = 12.9;   % INPUTE SOLENOID FIELD

%% Parameters 2
%mu0 = 4*pi*10^-7;
%Isol = B1*L/mu0; % Solenoid Current
%cJ= Isol/(Jmax*L); % HTS Thickness
%li = 0.67;   % Internal Inductance

%% Flux Equation
F(1) = mu0*R0*li*I/2 - pi*B1*x(1)^2 - pi*B1*( (x(2)^2)/6 + x(1)*x(2)/2 ); 

%% Stress Equation
F(2) = Sysol - pi*B1*Isol/(2*Lsol*x(2)*(x(2)-cJ))*(x(2)^2/6+x(2)*x(1)/2) - ...
    mu0*Isol^2/(4*Lsol^2*x(2)^2)*(x(2)^4/6 + 2/3*x(1)^3*x(2) + x(1)^2*x(2)^2)*...
    1/((x(1)+x(2))^2 - (0.5*(2*x(1)+x(2))+cJ/2)^2 + (0.5*(2*x(1)+x(2))-cJ/2)^2-x(1)^2);

end