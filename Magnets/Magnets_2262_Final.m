clear all;
close all;
clc;
clear;

global Sysol Lsol B1 Isol cJ mu0 R0 I li

%% General Input Parameters
a = 2;   % INPUT MINOR RADIUS
b = 1.19; % INPUT BLANKET THICKNESS
R0 = 6.3;  % INPUT MAJOR RADIUS
B0 = 5.2; % INPUT FIELD ON AXIS
Q_max = 350000; % INPUT HEAT FLUX [W/m^3]
N = 18; % NUMBER OF COILS
na = 1.5e20; % DENSITY OF PLASMA
Ta = 15e3*11600; % TEMPERATURE OF PLASMA
I = 15e6;    % INPUT PLASMA CURRENT

%% Parameters - Shape and TF Coils
Sy = 1050*10^6; % Maximum Allowable Stress TF
l0=0.722;      % Length factor
delta = 0.45;  % Triangularity
f=0.75;        % Fraction of straight magnet
kappa = 1.8;   % Elongation
mu0 = 4*pi*10^-7; % Permeability
eb = (a+b)/R0; % Normalized Inner Thickness
x0 = a+b;      % Minor Radius and Blanket Thickness
B_coil_max = B0/(1-eb);
ITF = 2*pi*B0*R0/mu0; % Total Current for given B0
Ic = ITF/18; % Current per coil
kB = 1.38*10^-23;
Price_HTS = 36; % Estimated cost of HTS per meter
Price_St = 9.6; % Estimated cost of steel per kg
Price_Cu = 8.3; % Estimated cost of copper per kg
Frac = 0.15; % Estimated HTS fraction in PF & CS
%% Solenoid input parameters
Jmax = (75)*10^6; % INPUT MAXIMUM ALLOWABLE CURRENT DENSITY PF
Jsol = 75e6; % LOWER LIMIT SET BY B1 - for 12.9 limit is 51
Sysol = 660e6;  % INPUT MAXIMUM ALLOWABLE STRESS
Lsol = 2*(kappa*a+b)*1.15;         % INPUT LENGTH OF SOLENOID
B1 = 12.9;   % INPUTE SOLENOID FIELD

%% TF Structure Calculations
% Parameterization and Cubic Solutions
k1 = (-(1+delta)+2*l0^4*(3-2*l0^2))/(l0^2*(1-l0^2)^2);
k2 = (2*(1+delta)-2*l0^2*(3-l0^4))/(l0^2*(1-l0^2)^2);
k3 = (-(1+delta)+2*l0^2*(2-l0^2))/(l0^2*(1-l0^2)^2);

c1 = (3/2)*(1/l0);
c2 = (1/2)*(1/l0^3);

% Cubic Solution
polynz = [1 k2/k3 k1/k3 (R0)/((a+b)*k3) + 1/k3];
z = roots(polynz);
z1 = z(1); z2 = z(2); z3 = z(3);

% Force Calculation

FR1 = -B0^2*R0^2/(2*mu0)*f*(kappa*a+b)/(R0-(a+b)); % Inward Centering

FR2 = -B0^2*R0^2/(2*mu0)*... % Outward Centering
    (c1/k3)*(kappa*a+b)/(a+b)*((1-3*c2/c1*z1)*acoth(sqrt(z1))/((z1-z2)*(z1-z3)*sqrt(z1))+...
     (3*c2/c1*z2-1)*acoth(sqrt(z2))/((z1-z2)*(z2-z3)*sqrt(z2))+...
     (3*c2/c1*z3-1)*acoth(sqrt(z3))/((z1-z3)*(z3-z2)*sqrt(z3)));

% Net Centering Force
 FC = 2*abs(FR1 + FR2);

% Tensile Force
 FZ = (pi*B0^2*R0^2/(mu0))*log((1+eb)/(1-eb));
 FT = (pi*B0^2*R0^2/(2*mu0))*log((1+eb)/(1-eb));

% Total Stress and Thickness
 c = (1/Sy)*(R0*FT/(pi*R0^2*(2-2*eb)) + FC/(2*f*(kappa*a + b))); % Thickness of material
 % cJ = R0*(1 - eb - ((1-eb)^2 - 2*B0/(mu0*R0*Jmax))^0.5);

% Arc Length Calculation
L1 = f*(a*kappa + b); % Straight Section Arc Length

func = @(x) sqrt((a + b).^2.*(2*k1*x + 4*k2*x.^3 + 6*k3*x.^5).^2 + ...
    (kappa.*a + b).^2.*(c1 - 3*c2*x.^2).^2);
L2 = integral(func,0,1); % Curved Section Arc Length

% Volume of TF Coil Structure not including HTS
V_TF = N*c*(2*pi/N)*(R0 - a - b - c/2)*(2*L1 + 2*L2);  % Volume of TF coil
C_TF = V_TF*Price_St*8000;

%% Winding Pack Calculations

Tape_w = 0.012; % tape width [m]
Tape_t = 0.0000446; % tape thickness [m]
Area_Tape = Tape_w*Tape_t;
Coil_P = 2*L1 + 2*L2; % Coil peremeter [m]
FOS = 0.6; % Stability margin

%Determine Tape Dimensions
Tape_N = ceil(Tape_w/Tape_t); % # of tapes per 12x12 mm HTS stack
Area_HTS = Tape_w^2; % Area of HTS in winding pack

%Determine Winding Pack Dimensions
It = 4670*B_coil_max^-0.796; % Critical current of 12mm tape @ 20K [A/mm2]
Jc = It/Tape_w/Tape_t/1000000; % Critical current density @ 20 K [A/mm2]
Jop = Jc*FOS; % Operating critical current density [A/mm2]
Imax = Jop*Tape_w^2*1000; % Max current per 12x12mm HTS stack [kA]
Turns1 = ceil(Ic/Imax/1000); % Number of cable turns per coil;
WP_AR = 2; % Winding pack aspect ratio
WP_d = ceil((Turns1/WP_AR)^0.5); % Winding pack radial depth [m]
WP_w = ceil(Turns1/WP_d); % Widing pack toroidal width [m]
Turns = WP_d*WP_w;
Cable_L = WP_d*Coil_P; % length of one cable [m]
Cable_N = WP_w; % # of cables per coil
Tape_L = Tape_N*Cable_L; % Length of tape per cable [m];

%Determine Copper Stablizer Dimensions
Tcs = 60; % lowest current sharing temperature
R1 = R0-a-b-c/2; % Radius from center of CS to TF straight leg
R2 = R0+a+b+c/2; % Radius from center of CS to TF outter D
Ro = (R1*R2)^0.5;
K = 0.5*log(R2/R1);
Coil_Inductance = (mu0*Turns^2*Ro*K^2)*(besseli(0,K) + 2*besseli(1,K) + besseli(2,K))/2; % TF coil inductance
Coil_Energy = 0.5*Coil_Inductance*(Imax*1000)^2; % TF coil magnetic energy
Cu_Cp = @(T) -0.000015917952141*T^4 + 0.001659188093121*T^3 - 0.003303875802231*T^2 - 0.076191481557999*T; % Specific heat of copper as a function of temperature
Cu_Rho = 8950; % Density of copper
Area_Cu = (Coil_Energy/Cu_Cp(Tcs)/Cu_Rho/(150-Tcs)/Coil_P/Turns)-(0.94*Area_HTS); % Area of Copper per cable

%Determine Steel Structure Dimensions
Force_St = Imax*Cable_L*B_coil_max; % estimate of force on cable due to magnetic field
St_Yield = 800 *10^6; % estimated yield strength of steel
Area_St = Imax*Cable_L*B_coil_max/((2/3)*St_Yield); % area of steel component of cable

%Determine Cooling Channel Dimensions
A = (Area_St+Area_Cu);
Cp_H2 = 10000; % Average specific heat of supercritical hydrogen @ 30bar from 20-55K
dT_H2 = 10; % Allowable temperature rise in the fluid
Q_int = (Q_max/0.029)*(1-exp(-.029*(WP_d))); % Integrated heat flux along one cable length [W/m^2]
Q_H2 = Q_int*(A); % Heat required to be removed from winding pack
Cable_mdot = Q_H2/Cp_H2/dT_H2; % Required mass flow rate to remove winding pack heat load per coil
f=0.015; % Cooling channel friction factor
dP_max = 0.01*10^6; % Allowable pressure drop in the cable
Cable_Dh = (f*Cable_L*(Cable_mdot)^2*8/(69.7*pi^2*(dP_max)))^(1/5); % finds cooling diameter based on mdot and max pressure drop allowed
Area_H2 = (pi/4)*Cable_Dh^2; % area of H2 cooling

%Determine Total Cable Cross Secional Area
Area_Cable = (Area_Cu + Area_HTS + Area_H2 + Area_St); % area of total cable
Cable_dim = Area_Cable^0.5; % square dimension of cable cross section

Frac_Cu = Area_Cu/Area_Cable; % fraction of copper
Frac_HTS = Area_HTS/Area_Cable; % fraction of HTS
Frac_H2 = Area_H2/Area_Cable; % fraction of H2
Frac_St = Area_St/Area_Cable; % fraction of steel

WP_Area = (WP_d*Cable_dim)*(WP_w*Cable_dim); % total winding pack area
Jeff = Imax/(Cable_dim^2*1000); % effective critical current density for entire winding pack

%Calculate Total WP Volumes and Costs
Vol_HTS = WP_Area*Frac_HTS*N*Coil_P;
Vol_Cu = WP_Area*Frac_Cu*N*Coil_P;
Vol_H2 = WP_Area*Frac_H2*N*Coil_P;
Vol_St = WP_Area*Frac_St*N*Coil_P;
Vol_WP = WP_Area*Coil_P*N;

Cost_Tape = Price_HTS*Tape_L*Cable_N;
Cost_Cu = Price_Cu*Vol_Cu*8950;
Cost_St = Price_St*Vol_St*8000;

%Calcualte Total TF Volume and Cost
Vol_TF = (Vol_WP)+V_TF;
Cost_TF = (Cost_Tape+Cost_Cu+Cost_St)*N+C_TF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POLOIDAL FIELD COIL Dimensions and Requirements
Bta = mu0*I/(2*pi*a);
% Pressure Calculation
nt=1.3;
nn = 0.5;
sigma=0.1;
xaa=0.3*(1-delta^2);
c0 = -delta/2;
c1 = 1/8*(9-2*delta-xaa);
c2 = delta/2;
c3 = 1/8*(-1+2*delta+xaa);
press = 2*na*Ta*kB*(1+nn)*(1+nt)/(1+nn+nt)*(c1 + a/R0*(6*c0*c1+3*c2*c3+c1...
   *(3+nn+nt)*(c2+sigma*(nn+nt)))/((2+nn+nt)*(3+nn+nt)));
betap = 2*mu0*press/(Bta^2);

% Calculate B Vertical
% Define Internal Inductance
rm = 0.6;
alpha = 1/(1-rm^2);

int1 = 1/4*(exp(2.*alpha).*(2.*alpha-1)+1);
int2 = (1+alpha).*(1-exp(2.*alpha));
int3 = 2*(1+alpha).*(exp(alpha)-1);
int4 = 27;

li = 1./(exp(alpha)-1-alpha).^2.*(int1+int2+int3+int4);
Bv = mu0*I/(4*pi*R0)*(betap + (li-3)/2 + log(8*R0/a));

 % Bz Field Contribution
 Bzfunc = @(rho,eta,B0,k) (B0/pi)*(1/((1+rho)^2+eta^2)^0.5)*(ellipticK(k)...
     + (1-rho^2-eta^2)/((1-rho)^2+eta^2)*ellipticE(k));

% Divertor Field Contribution
r1 = R0-a*(0.709*kappa+delta); z1 = 1.464*a*kappa; I1 = 0.485*I;
r2 = R0+a*(0.962*kappa-delta); z2 = a*kappa; I2 = -0.558*I;
r3 = R0+a*(0.962*kappa-delta); z3 = 1.272*a*kappa; I3 = 0.669*I;
r4 = R0-a*(0.709*kappa+delta); z4 = -1.464*a*kappa; I4 = 0.485*I;
r5 = R0+a*(0.962*kappa-delta); z5 = -a*kappa; I5 = -0.558*I;
r6 = R0+a*(0.962*kappa-delta); z6 = -1.272*a*kappa; I6 = 0.669*I;

rd = [r1 r6 r2 r5 r3 r4];
zd = [z1 z6 z2 z5 z3 z4];
Id = -[I1 I6 I2 I5 I3 I4];

B0d = (mu0/2)*Id./rd;
rhod = R0./rd;
etad = zd./rd;

kd=(4*rhod./((1+rhod).^2+etad.^2));

BzDiv = zeros(1,6);
for i = 1:6
BzDiv(i) = Bzfunc(rhod(i),etad(i),B0d(i),kd(i));
end

BzDivs = sum(BzDiv); % Total Divertor Field at Center of Plasma

% Poloidal Field Contribution
rp2 = R0 + a*(1.56*kappa - delta); zp2 = 1.35*a*kappa; Ip2 = -0.1750*I;
rp5 = R0 + a*(1.56*kappa - delta); zp5 = -1.35*a*kappa; Ip5 = -0.1750*I;
rp3 = R0 + a*(1.84*kappa - delta); zp3 = 0.83*a*kappa; Ip3 = -0.4250*I;
rp4 = R0 + a*(1.84*kappa - delta); zp4 = -0.83*a*kappa; Ip4 = -0.4250*I;

rp = [rp2 rp5 rp3 rp4];
zp = [zp2 zp5 zp3 zp4];
Ip = -[Ip2 Ip5 Ip3 Ip4];

B0p = (mu0/2)*Ip./rp;
rhop = R0./rp;
etap = zp./rp;

% m=k^2
kp = (4*rhop./((1 + rhop).^2 + etap.^2));

BzPF = zeros(1,4);
for i = 1:4
BzPF(i) = Bzfunc(rhop(i),etap(i),B0p(i),kp(i));
end

BzPFs = sum(BzPF); % Total PF Field at Center of Plasma

% Current Adjustment for Required Vertical B Field
Bvtotal = BzDivs + BzPFs;         % Total Field at center from Div and PF
scale = Bvtotal/Bv;               % With safety factor find ratio for new size reactor
Im=I;
ItPF = [Id Ip];

if scale < 1
 ItPF = 1.25*ItPF.*scale;   % Scale the currents with safety factor 1.25
 else if scale > 1
 ItPF = 1.25*ItPF./scale;
      end
end

% Calculate New BV due to PF and Divertor
BzDivnew = zeros(1,6);
B0dnew = (mu0/2)*ItPF(1:6)./rd;
for i = 1:6
    BzDivnew(i) = Bzfunc(rhod(i),etad(i),B0dnew(i),kd(i));
end
BzDivsnew = sum(BzDivnew);

BzPFnew = zeros(1,4);
B0pnew = (mu0/2)*ItPF(7:10)./rp;
for i = 1:4
    BzPFnew(i) = Bzfunc(rhop(i),etap(i),B0pnew(i),kp(i));
end
BzPFsnew = sum(BzPFnew);
Bvtotnew_SF = BzDivsnew + BzPFsnew;


% Forces on PF and Divertor Coils
% Axial Force

Fz = @(ri,rj,zi,zj,Ii,Ij) mu0*Ii*Ij/2*((zi-zj)/(sqrt((ri+rj)^2+(zi-zj)^2)))*...
    (ellipticK(sqrt(4*ri*rj/((ri+rj)^2+(zi-zj)^2)))...
    -ellipticE(sqrt(4*ri*rj/((ri+rj)^2+(zi-zj)^2)))*...
    (1+(1/2)*(4*ri*rj/((ri+rj)^2+(zi-zj)^2))/(1-4*ri*rj/((ri+rj)^2+(zi-zj)^2))));

% Radial Force
Fr = @(ri,rj,zi,zj,Ii,Ij) mu0*Ii*Ij/2*((ri+rj)/(sqrt((ri+rj)^2+(zi-zj)^2)))*...
    (ellipticK(sqrt(4*ri*rj/((ri+rj)^2+(zi-zj)^2)))*ri/(ri+rj)...
    + ellipticE(sqrt(4*ri*rj/((ri+rj)^2+(zi-zj)^2)))*1/(1-4*ri*rj/((ri+rj)^2+(zi-zj)^2)))*...
    ((4*ri*rj/((ri+rj)^2+(zi-zj)^2))/2 - ri/(ri+rj));

% Radial Self Force
Frs = @(Ii) mu0*(Ii^2)/2*(log(8*R0/a) - 1);

% Force due to plasma
Fzp = @(ri,rm,zi,zm,Ii,Im) mu0*Ii*Im/2*((zi-zm)/(sqrt((ri+rm)^2+(zi-zm)^2)))*...
    (ellipticK(sqrt(4*ri*rm/((ri+rm)^2 + (zi-zm)^2)))...
    -ellipticE(sqrt(4*ri*rm/((ri+rm)^2+(zi-zm)^2)))*...
    (1+(1/2)*(4*ri*rm/((ri+rm)^2+(zi-zm)^2))/(1-4*ri*rm/((ri+rm)^2+(zi-zm)^2))));

Frp = @(ri,rm,zi,zm,Ii,Im) mu0*Ii*Im/2*((ri+rm)/(sqrt((ri+rm)^2+(zi-zm)^2)))*...
    (ellipticK(sqrt(4*ri*rm/((ri+rm)^2+(zi-zm)^2)))*ri/(ri+rm)...
    + ellipticE(sqrt(4*ri*rm/((ri+rm)^2+(zi-zm)^2)))*1/(1-4*ri*rm/((ri+rm)^2+(zi-zm)^2)))*...
    ((4*ri*rm/((ri+rm)^2+(zi-zm)^2))/2 - ri/(ri+rm));

% Choose

Fzc = zeros(1,10);
Frc = zeros(1,10);

rm = 3.3;
zm = 0;

rp = [r1 r2 r3 r4 r5 r6 rp2 rp5 rp3 rp4 rm];
zp = [z1 z2 z3 z4 z5 z6 zp2 zp5 zp3 zp4 zm];
ItPF = [ItPF Im];
k = [kd kp];

for i = 1:10
   for j = 1:10
       if j == i
        Fzc(1,j) = Fzc(1,j)+0;
        Frc(1,j) = Frc(1,j)+Frs(ItPF(i));
       else
        Fzc(1,j) = Fzc(1,j) + Fz(rp(i),rp(j),zp(i),zp(j),ItPF(i),ItPF(j));
        Frc(1,j) = Frc(1,j) + Fr(rp(i),rp(j),zp(i),zp(j),ItPF(i),ItPF(j));
       end
   end
end

for i = 1:10
    Fzc(1,i) = Fzc(1,i) + Fzp(rp(i),rp(11),zp(i),zp(11),ItPF(i),ItPF(11));
    Frc(1,i) = Frc(1,i) + Frp(rp(i),rp(11),zp(i),zp(11),ItPF(i),ItPF(11));
end

% % Thickness of Steel for Vertical Poles of PF
% Fzc = abs(Fzc);
% Sy=600*10^6;
% rvz = sqrt((Fzc(7)+Fzc(8))/(2*pi*Sy));
%
% % Thickness of Radial Poles PF
% N=2;
% Jmax = 500*10^6;
% cmpoles = sqrt((abs(Fzc)/(N*pi*Sy)));

% Thickness of the PF Coil

Lp=0.3;
cj = zeros(1,10);
a2PF = zeros(1,10);
cm = zeros(1,10);
for p = 1:10
    cj(p) = abs(ItPF(p))/(Lp*Jmax);
    a2PF(p) = (abs(Frc(p))/2)/(Lp*Sy) + rp(p)+cj(p);
    cm(p) = a2PF(p)-rp(p) - cj(p);
end


%% Volume of PF Coils (HTS and Structure)

Vsc_PF = zeros(1,10);
Vst_PF = zeros(1,10);
for i = 1:10
Vsc_PF(i) = pi*(((rp(i) + a2PF(i))/2 + cj(i)/2).^2 - ((rp(i) + a2PF(i))/2 - cj(i)/2).^2)*Lp;
Vst_PF(i) = pi*((a2PF(i)).^2 - ((rp(i) + a2PF(i))/2 + cj(i)/2).^2)*Lp + ...
    pi*(((rp(i) + a2PF(i))/2 - cj(i)/2).^2 - (rp(i)).^2)*Lp;
end

Vol_PF = sum(Vsc_PF) + sum(Vst_PF);
Cost_PF = Frac*sum(Vsc_PF)*(Price_HTS/Area_Tape) + (((1-Frac)*sum(Vsc_PF))+sum(Vst_PF))*Price_St*8000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CENTRAL SOLENOID Dimensions and Requirements
%% Parameters

Isol = B1*Lsol/mu0; % Solenoid Current
cJ= Isol/(Jsol*Lsol); % HTS Thickness

%% Call fsolve to find roots of system of Flux & Stress Equations
func = @rootsol;
x0 = [0.4,0.2];
x = fsolve(func,x0,optimset('MaxFunEvals',100000000,'MaxIter',10000000));

%% Define Inner Radius and Total Thickness
a1 = x(1);
da = x(2);
%% Solenoid Thicknesses and Volumes
 cMCS = da - cJ; % Thickness of Structure
 a2CS = a1+da; % Outer Radius
 VJ_CS = pi*( ((a2CS+a1)/2 + cJ/2)^2 - ((a2CS+a1)/2 - cJ/2)^2)*Lsol; % Volume of HTS
 VM_CS = pi*( (a2CS)^2 - ((a2CS+a1)/2 + cJ/2)^2 + ((a2CS+a1)/2 - cJ/2)^2 - a1^2)*Lsol; % Volume of Structure


Vol_CS = VJ_CS + VM_CS;
Cost_CS = Price_St*(VM_CS+(1-Frac)*VJ_CS)*8000 + Frac*VJ_CS*(Price_HTS/Area_Tape);

Vol_ST_Total = (VM_CS + (sum(Vst_PF)) + (sum(Vsc_PF)*(1-Frac)) + VJ_CS*(1-Frac) + V_TF + Vol_WP*(1-Frac_HTS));
Cost_ST_Total = Vol_ST_Total*Price_St*8000;

Vol_HTS_Total = (VJ_CS*Frac+(sum(Vsc_PF)*Frac)+Vol_WP*Frac_HTS);
Cost_HTS_Total = Vol_HTS_Total*Price_HTS/(Area_Tape);

Cost_Total = Cost_HTS_Total+Cost_ST_Total;
Cost_Total11 = Cost_HTS_Total + Cost_ST_Total;

%end
% figure
% hold on
% plot(5:15,Cost_Total(5:15))
% plot(5:15,Cost_HTS_Total(5:15))
% plot(5:15,Cost_ST_Total(5:15))
% hold off
