%#######################Calculate the efficiency in Ehst and Karney 1991 paper
clear all;
close all;
clc;
%###########  Physical Constants
global e0 ee me mi mekev clight lnlambda Zeff R0 a epsilon B0 theta gamma0 nu_n nu_T n_bar T_bar coeff eta_LH n_para_m rho_J_m ;
e0 =  8.85e-12 ; %F/m
ee = 1.6021766208e-19 ; 
mekev = 510.998; %kev/c^2
me =  9.10938356e-31; %kg
mi = (3.02+2.01)/2.0*1.66e-27;  %D-T plasma, in kg
clight = 2.99792458e8; %light speed
%################### Plasma Parameters##############
lnlambda = 19.0; 
Zeff = 1.0;
R0 = 4.0; %major radius
a = 1.0; %minor radius
epsilon = a/R0;
B0 = 10; 
nu_n = 0.3; %density coefficient
nu_T = 1.2; %temperature cofficient
n_bar = 0.86; %n20
T_bar = 17.8;  %kev
theta = 0;   %the direciton to launch the wave
gamma0 = 8.562;
coeff = 1.0;
eta_LH = 0.75;
delta_n_para = 0.2;


%%%%%%%%%%%    Main part   %%%%%%%%%%%%%%%%%%%%%%%%%%%
solution = Solve();
n_para = solution(1);
rho_J = solution(2);
omega_nor2 = solution(3);
omega = solution(4);
exitflag = solution(5);
realflag = solution(6);

%%%%%%%%%%%%%%%%%%%%  Define rho_J_m ############
n_para_m = n_para+delta_n_para;
rho_J_m = sqrt( 1-(1-rho_J^2)*(1-delta_n_para/n_para_m )^(2.0/nu_T) );

efficiency = Efficiency([n_para_m, rho_J_m, omega]);
eta_CD = efficiency(1);
%eta_I  = efficiency(2);


%%%%%%%%%%%   Function to solve for n_para, rho_J, omega_nor2
function output = Solve()
global e0 ee me mi mekev clight lnlambda Zeff R0 a epsilon B0 theta gamma0 nu_n nu_T n_bar T_bar coeff eta_LH;

%coordinate information
xx = @(rho, theta)rho.*cos(theta);
R = @(rho,theta)R0+a*rho.*cos(theta);
B = @(rho,theta)B0./(1+epsilon*rho*cos(theta));
BM = @(rho)B(abs(rho),pi); 

n_profile = @(rho)n_bar*(1+nu_n)*(1-rho.^2).^nu_n; %density profile in radial flux label
T_profile = @(rho)T_bar*(1+nu_T)*(1-rho.^2).^nu_T; %temperature profile in radial flux label
n_max = n_bar*(1+nu_n);
T_max = T_bar*(1+nu_T);
v_Te = @(rho)sqrt(2*T_profile(rho)./mekev)*clight;
v_e = @(rho)sqrt(T_profile(rho)./mekev)*clight;

omega_pe2 = @(rho,theta)(n_profile(rho)*1e20).*ee^2./(me*e0); %plasma frequency square

omega_ce = @(rho,theta)ee*B(rho,theta)/me; %electron cyclotron frequency
omega_ci = @(rho,theta)ee*B(rho,theta)/mi; %ion cyclotron frequency

X = @(rho,theta)omega_pe2(rho,theta)./omega_ce(rho,theta).^2;  %(eqn5)

%x_initial = [2.6569,0.0505,3.2456];
x_initial = [1,0.5,0.5];
[sol,fval,exitflag] = fsolve(@Eqnset,x_initial);

if exitflag >=1
	exitflag = 1;
else 
	exitflag = 0;
end

n_para = sol(1);
rho_J = sol(2);
omega_nor2 = sol(3);
omega = sqrt(omega_nor2.*omega_ce(rho_J,theta).*omega_ci(rho_J,theta));   %(eqn15)
%%% check if it is real solution
if abs(imag(n_para))<1e-7 & abs(imag(rho_J))<1e-7 & abs(imag(omega_nor2))<1e-7
	realflag = 1;
else
	realflag = 0;
end

output(1) = n_para;
output(2) = rho_J;
output(3) = omega_nor2;
output(4) = omega;
output(5) = exitflag;
output(6) = realflag;

end

%%%%%%%%%%%    Equations to solve %%%%%%%%%%%%%%%%%%%55
function F = Eqnset(vars)
global e0 ee me mi mekev clight lnlambda Zeff R0 a epsilon B0 theta gamma0 nu_n nu_T n_bar T_bar coeff eta_LH;
xx = @(rho, theta)rho.*cos(theta);
R = @(rho,theta)R0+a*rho.*cos(theta);
B = @(rho,theta)B0*R0./R(rho,theta);
BM = @(rho)B(abs(rho),pi); 
n_profile = @(rho)n_bar*(1+nu_n)*(1-rho.^2).^nu_n; %density profile in radial flux label
T_profile = @(rho)T_bar*(1+nu_T)*(1-rho.^2).^nu_T; %temperature profile in radial flux label
omega_pe2 = @(rho,theta)(n_profile(rho)*1e20).*ee^2./(me*e0);
omega_ce = @(rho,theta)ee*B(rho,theta)/me;
omega_ci = @(rho,theta)ee*B(rho,theta)/mi;
%########### Wave Parameters###########
X = @(rho,theta)omega_pe2(rho,theta)./omega_ce(rho,theta).^2;
%%%%%%%%%%%%%%%%   (eqn12)
n_para = vars(1);
rho = vars(2);
omega_nor2 = vars(3);
F(1) =  coeff*n_para.^2 - ( (1-(1-omega_nor2).*X(rho,theta)./omega_nor2 ).^0.5 + X(rho,theta).^0.5 ).^2;
F(2) = omega_nor2 - 0.5*X(rho,theta)./(1+X(rho,theta)) - 0.5*( X(rho,theta).^2./(1+X(rho,theta)).^2 + 4*gamma0^2*X(rho,theta).^3./(1+X(rho,theta)) ).^0.5;
F(3) = (1+nu_T)*(1-rho.^2).^nu_T.*n_para.^2 - 28.4./T_bar;
end

%%%%%%%%%%%%   Function to calculate the efficiency
function output = Efficiency(vars)
global e0 ee me mi mekev clight lnlambda Zeff R0 a epsilon B0 theta gamma0 nu_n nu_T n_bar T_bar coeff eta_LH n_para_m ;
xx = @(rho, theta)rho.*cos(theta);
R = @(rho,theta)R0+a*rho.*cos(theta);
B = @(rho,theta)B0*R0./R(rho,theta);
BM = @(rho)B(abs(rho),pi); 
n_profile = @(rho)n_bar*(1+nu_n)*(1-rho.^2).^nu_n; %density profile in radial flux label
T_profile = @(rho)T_bar*(1+nu_T)*(1-rho.^2).^nu_T; %temperature profile in radial flux label
v_Te = @(rho)sqrt(2*T_profile(rho)./mekev)*clight;
v_e = @(rho)sqrt(T_profile(rho)./mekev)*clight;
%########### Wave Parameters###########
X = @(rho,theta)omega_pe2(rho,theta)./omega_ce(rho,theta).^2;
n_para_m = vars(1);
rho_J_m = vars(2);
omega = vars(3);

%%%%%%%%%%%%%%%   Calculating efficiency   %%%%%%%%%%%%%%
ww = clight./(v_e(rho_J_m).*n_para_m);   %(eqn6)
% following are (eqn5)
xt2 = ww.^2*(B(rho_J_m,theta)./(BM(rho_J_m)-B(rho_J_m,theta)));
MM = 1;
nn = 0.77;
xr = 3.5;
RR = 1- epsilon^nn*rho_J_m^nn*(xr^2+ww^2)^0.5/(epsilon^nn*rho_J_m^nn*xr+ww);
mm = 1.38;
cc = 0.389;
CC = 1-exp(-cc^mm*xt2^mm);
KK = 3.0/Zeff;  
DD = 3.83/Zeff^0.707;
%%%%%%   efficiency in different units  %%%%%%%%%%%%%
eta0 = KK/ww+DD+4*ww^2/(5+Zeff);
eta_tilde = CC*MM*RR*eta0; %(eqn13)
eta = (0.384/lnlambda)*T_profile(rho_J_m)/n_profile(rho_J_m)*eta_tilde;  %(eqn2)
eta_I = eta/(2*pi*R(rho_J_m,theta));
% eta_CDv2 = eta_I*n_profile(rho_J_m)*R(rho_J_m,theta);
eta_CD = (0.06108/lnlambda)*eta_LH*T_profile(rho_J_m)*eta_tilde;  %(eqn14)
output(1) = eta_CD;
%output(2) = eta_I;
end

