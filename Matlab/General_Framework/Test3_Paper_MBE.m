% TEST 3 of the Paper: AOI WITH HETEROGENEOUS MIXING (part of EPITIME software)
%
% E. Messina, M. Pezzella, A. Vecchio.
% "Nonlocal finite difference discretization of a class of renewal equation
%  models for epidemics"
% Mathematical Biosciences and Engineering, 2023, 20(7): 11656-11675.
% DOI: 10.3934/mbe.2023518
%
% For further details on the Age of Infection model with heterogeneous  
% mixing we refer also to:
% E. Messina, M. Pezzella, A. Vecchio.
% "A long-time behavior preserving numerical scheme for age-of-infection
%  epidemic models with heterogeneous mixing"
% Applied Numerical Mathematics, 2024, 200: 344-357.
% DOI: 10.1016/j.apnum.2023.04.009
%
% ==========================================================================
% Version:    1.0
% Date:       April, 2026
% Authors:    M. Pezzella, G. Zanghirati
% See:        https://arxiv.org/abs/2605.00067v1
% License:    GNU GPL v.3
% Repository: https://github.com/ghitan/EPITIME


%====================|
% Problem definition |
%====================|
%%% Two alternative input modes
%%% (1) Structured input:
%%% => NSFD_Renewal( Ren_Problem )
% N = [100; 900];
% Ren_Problem = struct();
% Ren_Problem.S0    = [99; 899] ; 
% Ren_Problem.A     = @(i,t)( normpdf(t, 0.2, 0.6) );
% Ren_Problem.phi0  = @(i,t)( (N(i) - S0(i)) * A(i,t) );
% Ren_Problem.P0    = @(t)( zeros(size(t)) );
% Ren_Problem.B     = @(t)( zeros(size(t)) ); 
% Ren_Problem.T     = 10;
% Ren_Problem.h     = 1e-1;
% Ren_Problem.beta  = [5; 10];
% Ren_Problem.BETA  = [0.4/N(1), 0.6/N(2); 0.5/N(1), 0.5/N(2)];
% Ren_Problem.alpha = zeros(2,1);
% Ren_Problem.c     = zeros(2,1);

% (2) Positional input:
% => NSFD_Renewal(S0, phi0, P0, A, B, beta, BETA, alpha, c, h, T)
N = [100; 900];        S0 = [99; 899];   beta = [5; 10];
A = @(i,t)( normpdf(t, 0.2, 0.6) );
h = 1e-1;              T = 10;           phi0 = @(i,t)( (N(i) - S0(i)) * A(i,t) );
P0 = @(t)( zeros(size(t)) );                B = @(t)( zeros(size(t)) );    
BETA = [0.4/N(1), 0.6/N(2); 0.5/N(1), 0.5/N(2)];
alpha = zeros(2,1);       c = zeros(2,1);           verbosity = 0;

[t,S,phi,P] = NSFD_Renewal(S0,phi0,P0,A,B,beta,BETA,alpha,c,h,T,verbosity);

%========================|
% Continuous  Final Size |
%========================|
Aint = arrayfun(@(i) integral(@(t) A(i, t), 0, inf), 1:2);
FS_Function = @(x) [ ...
    log( S0(1) / x(1) ) - 5  * ( 0.4 / N(1) * (N(1) - x(1)) * Aint(1) ...
        + 0.6 / N(2) * (N(2) - x(2)) * Aint(2) ); ...
    log( S0(2) / x(2) ) - 10 * ( 0.5 / N(1) * (N(1) - x(1)) * Aint(1) ...
        + 0.5 / N(2) * (N(2) - x(2)) * Aint(2) ) ];
options = optimoptions('fsolve', 'FunctionTolerance', 1e-14, ...
    'StepTolerance', 1e-14, 'OptimalityTolerance', 1e-14, ...  
    'Display', 'off');      
S_infty = fsolve(FS_Function, [1; 1], options);

%================================|
% Plots of the Numerical Soution |
%================================|
figure(1);
subplot(1,2,1)
plot(t,S(1,:),'k','Linewidth',2); hold on; 
plot(t,phi(1,:),'k--','Linewidth',2); axis square
plot(t(end),S_infty(1),'k.','MarkerSize',25);
legend('S₁(t)','φ₁(t)','S₁(∞)'); 
title(sprintf('Group 1, h = %.2e', h));       xlabel('t');
subplot(1,2,2)
plot(t,S(2,:),'k','Linewidth',2); hold on; 
plot(t,phi(2,:),'k--','Linewidth',2); axis square
plot(t(end),S_infty(2),'k.','MarkerSize',25);
legend('S₂(t)','φ₂(t)','S₂(∞)'); 
title(sprintf('Group 2, h = %.2e', h));       xlabel('t');      axis tight;

%============================================|
% Num. Final Size for different values of h  |
%============================================|
T = 100;  
N_test = 4; final_size_num = zeros(size(S,1), N_test);
error_norm = zeros(size(S,1), N_test);
for j = 1:N_test
    h = 2^(-j)/10;
    [~,S,~,~] = NSFD_Renewal(S0,phi0,P0,A,B,beta,BETA,alpha,c,h,T,verbosity);
    final_size_num(:,j) = S(:,end);
    tmp = S(:,end)-S_infty;
    error_S1 = S(1,end)-S_infty(1);
    error_S2 = S(2,end)-S_infty(2);
    error_norm(j) = norm(tmp);
end
order_norm = log2(error_norm(1:end-1)./error_norm(2:end));

%=================================================|
% Table with Num. Final Size and Absolute Errors  |
%=================================================|
fprintf('\n %-11s %-17s %-19s %-15s %-11s\n', ...
        'h', 'Abs Err. S₁(∞)', 'Abs Err. S₂(∞)', '‖Err(S)‖₂', 'Order');
for j = 1:4
    h_label = sprintf('0.1*2^-%d', j);
    err_S1 = final_size_num(1,j) - S_infty(1);
    err_S2 = final_size_num(2,j) - S_infty(2);
    err_norm_j = error_norm(j);
    if j == 1
        fprintf('%-12s %-18.5e %-18.5e %-15.5e %-10s\n', ...
                h_label, err_S1, err_S2, err_norm_j, '\\');
    else
        fprintf('%-12s %-18.5e %-18.5e %-15.5e %-10.4f\n', ...
                h_label, err_S1, err_S2, err_norm_j, order_norm(j-1));
    end
end




