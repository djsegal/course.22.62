clear 

%***************************************************************
% Calculate the bootstrap fraction 
%***************************************************************
% inputs 
a = 1;                     % minor radius 
R0 = 5;                    % major radiua 
mu_0 = 8.85e-12;           % permativity of free space (m^-3 kg^-1 s^4 A^2)
kappa = 1.76;              % triangularity 
a_bar = kappa^(1/2) * a;       % ellipse adjusted 


% assume density and temperature profiles 
nu_n = 0.3;     % density profile constant 
nu_T = 1.2;     % temperature profile constant 
n_bar = 1;      % avg. density 
T_bar = 5;      % avg. temperature 

n = @(rho) n_bar .* (1 + nu_n) .* (1 - rho.^2).^nu_n;
T = @(rho) T_bar .* (1 + nu_T) .* (1 - rho.^2).^nu_T; 

% gradient of temperature and density 
% used in bootstrap current formula (eqn 11 in Fusion Design 2) 
grad_n = @(rho) 2 .* n_bar .* (1 + nu_n) .* nu_n .* rho .* (1 - rho.^2).^(nu_n - 1);
grad_T = @(rho) 2 .* T_bar .* (1 + nu_T) .* nu_T .* rho .* (1 - rho.^2).^(nu_T - 1);

% Assume functional form of j_phi (total) 
alpha = 1/(1 - 0.6^2);
j_phi = @(rho) (alpha^2 * (1 - rho.^2) .* exp(alpha .* rho.^2))/...
                (exp(alpha) - alpha - 1 );            
b_theta = @(rho) (2./rho) .* integral(@(rho)j_phi(rho).*rho, 0,1); 


% Function to calculate Bootstrap Current (eqn 11 in Fusion Design 2) 
% B_theta left as parameter 
JB = @(rho)...
    4.88 .* (rho./R0).^(1/2) .* ((n(rho) .* T(rho))./b_theta(rho)) .* ...
    ((1./n(rho)) .* grad_n(rho) + 0.055 .* (1./T(rho)) .* grad_T(rho)); 


I_B = 2 * pi * a_bar^2 * integral(@(rho)JB(rho).*rho, 0, 1); 


