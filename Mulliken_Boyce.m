clc
close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------Mullken & Bryce constituative model for 1D uniaxial loading--------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auther : Asger Riis Vienberg
% 
% This matlab program impliments the model proposed by Mullken & Bryce
% DOI: 10.1016/j.ijsolstr.2005.04.016



%--------------------------------------------------------------------------
%   Inputs
%--------------------------------------------------------------------------

% Material properties (read off graf for theta = 300K)

    % elastisk properties
mat = struct();                                 
mat.E_alpha = 1700e6;
mat.E_beta = 700e6;
mat.nu = 0.38;
mat.mu_alpha = mat.E_alpha/(2*(1+mat.nu));
mat.mu_beta = mat.E_beta/(2*(mat.nu+1));
mat.K_alpha = mat.E_alpha/(3*(1-2*mat.nu));
mat.K_beta = mat.E_beta/3*(1-2*mat.nu);

    % Internal state varible s
s = struct;
s.alpha = 0.077*mat.mu_alpha/(1-mat.nu);
s.beta = 0.077*mat.mu_beta/(1-mat.nu);

    % nonlinear langvin spring
mat.C_R = 14.2e6;
mat.sqrtN = 2.3;
    
    % Plastic strain evloution constants for alpha component
mat.dot_gamma_0_alpha = 2.94e16;
mat.DeltaG_alpha = 3.744e-19;
mat.alpha_p_alpha = 0.168;
mat.h_alpha = 250e6;
mat.s_ss_alpha = 0.67*s.alpha;

    % Plastic strain evloution constants for beta component
mat.dot_gamma_0_beta = 3.39e5;
mat.DeltaG_beta = 3.769e-20;
mat.alpha_p_beta = 0.245;
mat.h_beta = 0;
mat.s_ss_beta = 0;

    % Enviromental and physical constants
p_ambient  = 1e5            % pressure, assumed to be the ambient
theta = 150;                % Temperture [K]
k = 1.380649e-23;           % Boltsmanns constant


% Solver settings
n = 20000;                  % Number of timesteps
    
% Testing parameters
totalTrueStrain =  0.8;     % total strain of the test
trueStrainRate = 5050;      % Strain rate duing the test


%--------------------------------------------------------------------------
%   end of inputs
%--------------------------------------------------------------------------


% Derived solver parameters
t_tot = totalTrueStrain/abs(trueStrainRate);    % Total time of test
dt = t_tot/n;                                   % Size of timesteps
t = 0;                                          % Initial time

% Intilize varibles

    % Total deformation gradient
F = diag([1 1 1]);                   

    % Plastic part of deformation
F_p = struct;                        
F_p.alpha = diag([1,1,1]);
F_p.beta = diag([1,1,1]);
    
    % strech rate
D = struct;
D.alpha = 0;
D.beta = 0;

    % plastic strain rate
dot_F_p = struct;
dot_F_p.alpha = 0;
dot_F_p.beta = 0;

    % equivalent effective shear stress
tau = struct;                                   
tau.alpha = 0;
tau.beta = 0;

    % magnitude of plastic strainrate
dot_gamma = struct;
dot_gamma.alpha = 0;
dot_gamma.beta = 0;
    
    % Direction of plastic strain rate
N = struct;
N.alpha = 0;
N.beta = 0;
    
    % evolution rate of internal state varible
dot_s = struct;
dot_s.alpha = 0;
dot_s.beta = 0; 


    % History varibles
strech_hist = ones(1,n+1);          % Strain history 
sigma_hist = zeros(n+1,1);          % Stress history of sigma_11
e_true = zeros(n+1,1);              % True strain history

% Time intergration 
for i = 1:n
    
    % Update time
    t = t + dt;
    trueStrain = trueStrainRate*t;

    % Update strain in testing direction
    F(1,1) = exp(trueStrain);

    % Update the deforemation gradient and find stress in the
    % satisfing that the total normal stress are = [sigma11 0 0 ].

    [T, F] = MB_solve(F, F_p, mat);

    % Evaluate state parameters and calculate the new evolution values

        % effective equvlent shear stress
    tau.alpha = sqrt(1/2*sum((T.alpha) .* (T.alpha),'all'));
    tau.beta = sqrt(1/2*sum((T.beta) .* (T.beta),'all'));

        % direction of plastic strain
    N.alpha = dev(T.alpha)/sqrt(sum((T.alpha) .* (T.alpha),'all'));
    N.beta = dev(T.beta)/sqrt(sum((T.beta) .* (T.beta),'all'));

        % Plastic strain rate
    p = 1/3*trace(T.total) + p_ambient;
    dot_gamma.alpha = mat.dot_gamma_0_alpha* exp(-mat.DeltaG_alpha/(k*theta) * (1 - tau.alpha/(s.alpha + mat.alpha_p_alpha*p_ambient)) );
    dot_gamma.beta = mat.dot_gamma_0_beta* exp(-mat.DeltaG_beta/(k*theta) * (1 - tau.beta/(s.beta + mat.alpha_p_beta*p_ambient)) );

        % evolution rate of internal varible
    dot_s.alpha = mat.h_alpha*(1 - s.alpha/mat.s_ss_alpha) * dot_gamma.alpha;
    dot_s.beta = 0; % mat.h_beta*(1 - s.beta/mat.s_ss_beta) * dot_gamma.beta;

        % plastik strech rate
    D.alpha = dot_gamma.alpha*N.alpha;
    D.beta = dot_gamma.beta*N.beta;

        % plastik stain rate
    dot_F_p.alpha = D.alpha*F_p.alpha;
    dot_F_p.beta = D.beta*F_p.beta;


    % Intergrating varibles in time

        % Plast strain
    F_p.alpha = F_p.alpha + dt*dot_F_p.alpha;
    F_p.beta = F_p.beta + dt*dot_F_p.beta;

        % Internal state varible
    s.alpha =s.alpha + dt*dot_s.alpha;
    s.beta =s.beta + dt*dot_s.beta;

    % Update history varibles
    e_true(i+1) = log(F(1,1));
    strech_hist(i+1) = F(1,1);
    sigma_hist(i+1) = T.total(1,1);
end


% ploting
plot(e_true, abs(sigma_hist)/1e6, '.-')
hold on
grid on

% Plot labels
xlabel('Stain')
ylabel('Stress[MPa]')




%% functions
function [T,F] = MB_solve(F,F_p,mat)

%under relaxation
relaxationFactor = 1;
maxIterations = 1000;

for k = 1:maxIterations

    %evaluat function value
    T = MB_stresses(F,F_p,mat);
    
    %pertubated varibles
    h = 1e-8;

    dF = F + h*diag([0,1i,1i]);

    %evaluate the complex step
    complexT = MB_stresses(dF, F_p, mat);

    %evaluat the derivative
    dT = imag(complexT.total)/h;

    %calculate the step size
    step = relaxationFactor * T.total(2,2)/dT(2,2);

    % update varibles

    F(2,2) = F(2,2) - step;
    F(3,3) = F(2,2);

    % evaluate convergens
    if abs(step) < 1e-12
        % calculate elastic deformations
        break
    end
end

T = MB_stresses(F,F_p,mat);

if k == maxIterations
    error('Stress solution not converged')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = MB_stresses(F, F_p, mat)

% calculate elastic deformations
F_e.alpha = F*inv(F_p.alpha);
F_e.beta = F*inv(F_p.beta);

lnU.alpha = logm(F_e.alpha);
lnU.beta = logm(F_e.beta);
%computes the 3 compontets of the MB stress
T = struct;

%computes the alpha componentes
T.alpha = mat.E_alpha*(lnU.alpha + mat.nu*trace(lnU.alpha)*eye(3)/(1 - 2*mat.nu))/(1 + mat.nu);
T.alpha = 1/det(F_e.alpha)*T.alpha;

%computes the beta component
T.beta = mat.E_beta*(lnU.beta + mat.nu*trace(lnU.beta)*eye(3)/(1 - 2*mat.nu))/(1 + mat.nu);
T.beta = 1/det(F_e.beta)*T.beta;

% computes the B componentes
B_bar = det(F)^(-2/3)*(F*F');
lambda_p_chain = sqrt(trace(B_bar)/3);

T.B = mat.C_R/3*mat.sqrtN*inv_L(lambda_p_chain/mat.sqrtN)*dev(B_bar)/lambda_p_chain;

T.total = T.alpha + T.beta + T.B;

end



function beta = inv_L(alpha)

beta = 0;
for i = 1:100

    if beta == 0
        L = 0 - alpha;
        dL = 1/3;
    else
        L = coth(beta) - 1/beta - alpha;
        dL = 1 - coth(beta)^2 + 1/beta^2;
    end

    beta = beta - L/dL;

    if abs(L/dL) < 1e-12
        break
    end

end

if i == 100
    error('Error: Invers Langevin did not converge')
end

end


%%%%%%%%%%%%%%%%%%%%%%%

function A = dev(B)
A = B - 1/3*trace(B)*eye(3);
end




