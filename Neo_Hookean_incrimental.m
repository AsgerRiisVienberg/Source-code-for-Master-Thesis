clc
clear all
close all

% This code solve the load case of uniaxial tension of a neo-hookean
% material using Euler interation. It plots it for a series of increasing
% loadsteps.
%
% For Euler intergration the update rule is 
% dsigma = Kt * dF
%
% Where dsigma is the change in stress for each loadstep, dF is the
% corrosponding change in the deformation gradient and Kt is the tangent
% stiffness matrix. Kt is defined as 
%
%  C_ijkl = d sigma_ij / d F_kl
% 
% For the case of Uniaxial tension this reduces to a 3x3 system of
% equations.
%
% Kt = d sigma_ij / d F_kl * korndelta_ij 

j = 1;
for n = [round(logspace(1,4,5))] % defind the number of loadsteps
    
    P_tot = 1.586100888e7;    % Total  load
    dP = P_tot/n;   % Size of loadstep
    P = 0;          % Initial load
    A = 1;          % Crossection area, set to unity for convience

    %Hyper-elastic constants in SI units 
    mu= 2e6;    
    K = 99e6;

    % Initilization of the deformation gradinet tangent stiffness
    F = diag(ones(1,3));    % sretches are 1 when stressfree
    Kt = zeros(3);
    
    % Initilize plotting varibles
    sigma_hist = zeros(n+1,1);  % Stress history
    lambda_hist = ones(n+1,1);  % Strecht history

    % Incrimentation loop
    for i = 1:n
        P = P + dP;             % Incriment load
        sigma = P/A;            % Calculate current stresslevel
        dsigma = [dP/A; 0; 0];  % Define increment in load vector
        
        % Evaluate stiffness matrix (see maple dokument for details
        Kt(1,1) = (0.222222222*mu*F(1,1)^2 + 0.5555555555*mu*F(2,2)^2 + 0.5555555555*mu*F(3,3)^2 + K*F(1,1)^2*F(2,2)^2*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3))/(F(2,2)*F(3,3)*F(1,1)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3));
        Kt(1,2) = (-1.111111111*mu*F(1,1)^2 - 0.1111111111*mu*F(2,2)^2 + 0.5555555555*mu*F(3,3)^2 + K*F(1,1)^2*F(2,2)^2*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3))/(F(1,1)*F(3,3)*F(2,2)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3));
        Kt(1,3) = (-1.111111111*mu*F(1,1)^2 - 0.1111111111*mu*F(2,2)^2 + 0.5555555555*mu*F(3,3)^2 + K*F(1,1)^2*F(2,2)^2*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3))/(F(1,1)*F(3,3)*F(2,2)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3));
        Kt(2,1) = (-1.111111111*mu*F(2,2)^2 - 0.1111111111*mu*F(1,1)^2 + 0.5555555555*mu*F(3,3)^2 + K*F(1,1)^2*F(2,2)^2*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3))/(F(2,2)*F(3,3)*F(1,1)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3));
        Kt(2,2) = (0.222222222*mu*F(2,2)^2 + 0.5555555555*mu*F(1,1)^2 + 0.5555555555*mu*F(3,3)^2 + K*F(1,1)^2*F(2,2)^2*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3))/(F(1,1)*F(3,3)*F(2,2)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3));
        Kt(2,3) = (-1.111111111*mu*F(2,2)^2 + 0.5555555555*mu*F(1,1)^2 - 0.1111111111*mu*F(3,3)^2 + K*F(1,1)^2*F(2,2)^2*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3))/(F(1,1)*F(2,2)*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3));
        Kt(3,1) = (-1.111111111*mu*F(3,3)^2 - 0.1111111111*mu*F(1,1)^2 + 0.5555555555*mu*F(2,2)^2 + K*F(1,1)^2*F(2,2)^2*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3))/(F(2,2)*F(3,3)*F(1,1)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3));
        Kt(3,2) = (-1.111111111*mu*F(3,3)^2 + 0.5555555555*mu*F(1,1)^2 - 0.1111111111*mu*F(2,2)^2 + K*F(1,1)^2*F(2,2)^2*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3))/(F(1,1)*F(3,3)*F(2,2)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3));
        Kt(3,3) = (0.222222222*mu*F(3,3)^2 + 0.5555555555*mu*F(1,1)^2 + 0.5555555555*mu*F(2,2)^2 + K*F(1,1)^2*F(2,2)^2*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3))/(F(1,1)*F(2,2)*F(3,3)^2*(F(1,1)*F(2,2)*F(3,3))^(2/3));

        dF = Kt\dsigma;     % Solve system

        F = F + diag(dF);   % Update deformation gradient

        % Save results to history varibles
        sigma_hist(i+1) = sigma;
        lambda_hist(i+1) = F(1,1);
    end
    lambda_err(j) = (lambda_hist(end) - 3)/3*100;
    legendList{j} = sprintf('%d',n);
    j = j + 1;
    % plot current solution
    figure(1)
    hold on
    plot(lambda_hist,sigma_hist,'.-')
    % axis([1 5 0 35])
    grid on
    


end
figure(1)
ylabel('Stress [MPa]')
xlabel('Strech \lambda_{11}')
legend(legendList,'Location','northwest')

leg = legend('show');
title(leg,'Number of Euler points')
enhance_plot()
hold off

figure(2)
loglog([round(logspace(1,4,5))],abs(lambda_err),'.-')
xlabel('Number of Euler points')
ylabel('Error [%]')
enhance_plot(0,0,0,0,-1)
grid on





