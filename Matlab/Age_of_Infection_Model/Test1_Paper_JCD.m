% TEST 1 of the Paper: EXPERIMENTAL ORDER OF CONVERGENCE (part of EPITIME software)
% 
% E. Messina, M. Pezzella, A. Vecchio.
% "A non-standard numerical scheme for an age-of-infection epidemic model"
% Journal of Computational Dynamics, 2022, 9(2): 239–252.
% DOI: 10.3934/jcd.2021029
% 
% ==========================================================================
% Version:    1.0
% Date:       April, 2026
% Authors:    M. Pezzella, G. Zanghirati
% See:        https://arxiv.org/abs/2605.00067v1
% License:    GNU GPL v.3
% Repository: https://github.com/ghitan/EPITIME

clearvars; clc; close all;
fprintf('\nEPITIME AoI: test 1: experimental order of convergence\n');

%====================|
% Problem definition |
%====================|
% Two alternative input modes
% (1) Structured input:
% => NSFD_AoI( AoI_problem) 
AoI_problem = struct();
AoI_problem.T = 1;
AoI_problem.h = 2^(-16);
AoI_problem.A = @(t)( (1+t).^(-2) ); 
AoI_problem.N = 1e1;
AoI_problem.beta = 1e-3;
AoI_problem.S0 = 0.9*AoI_problem.N;
AoI_problem.phi0 = @(t)( (AoI_problem.N - AoI_problem.S0) * AoI_problem.A(t) );
h_ref = AoI_problem.h;
AoI_problem.verbosity = 0;
%%% (2) Positional input:
%%% => NSFD_AoI(S0, N, T, h, beta, A, phi0, verbosity)
% T = 1;    N = 1e1;    S0 = 0.9*N;     beta = 0.3;
% A = @(t)( (1+t).^(-2) ); phi0 = @(t)( (N - S0) .* A(t) );
% h_ref = 2^(-16);    verbosity = 0;

%====================|
% Reference Solution |
%====================|
fprintf('\nComputing the reference solution...')
[~, Y_ref, P_ref] = NSFD_AoI( AoI_problem );
% [~, Y_ref, P_ref] = NSFD_AoI(S0, N, T, h_ref, beta, A, phi0);    % Alt. (2)
fprintf(['done.\nThe reference solution was computed in %.3f seconds with a' ...
         'stepsize h = 2^%d = %.2e.\n'], P_ref.elapsed_time, log2(h_ref),h_ref);

%================================|
% Experimental Convergence Study |
%================================|
start_expn = 6; TestN = -log2(h_ref) -2 -start_expn;
H = zeros(TestN,1);     S_err = zeros(TestN,1);   phi_err = zeros(TestN,1); 

% Prepare the running echo to console
Ndisp = max(1., TestN / 10.) * ( 1 : min(TestN, 10) )';
k = 1;
fprintf('\nConvergence study started...\n0%%..');

for j = 1:TestN
    h_test = 2^-(j + start_expn);
    H(j) = h_test;
    AoI_problem_test = AoI_problem;            AoI_problem_test.h = h_test;
    [~, Y_test, ~] = NSFD_AoI( AoI_problem_test );
    % [~,Y_test,~] = NSFD_AoI(S0, N, T, h_test, beta, A, Phi0);  % Alt. (2)
    zooming = 1:h_test/h_ref:size(Y_ref,2);
    S_err(j) = norm(Y_ref(1,zooming)-Y_test(1,:)) / norm(Y_test(1,:));
    phi_err(j) = norm(Y_ref(2,zooming)-Y_test(2,:)) / norm(Y_test(2,:));
    % Display running iterations percentage
    if ( j >= Ndisp(k) ), fprintf('%d%%..', fix(100. * Ndisp(k)/ TestN)); k = k + 1; end
end
fprintf('%c%c\nConvergence study finished.\n', char(8), char(8));

OrderS   = log2(S_err(1:end-1) ./ S_err(2:end));
OrderPhi = log2(phi_err(1:end-1) ./ phi_err(2:end));
OrderS   = [0; OrderS];
OrderPhi = [0; OrderPhi];

%==============================================|
% Relative Errors and Experimental Order Table |
%==============================================|
fprintf('\n%-10s %-15s %-10s %-15s %-10s\n', 'h', 'RelError_S', 'Order_S',...
        'RelError_phi', 'Order_phi');
for j = 1:length(H)
    fprintf('2^-%-7d %-15.3e %-10.2f %-15.3e %-10.2f\n', ...
            j + start_expn, S_err(j), OrderS(j), phi_err(j), OrderPhi(j));
end

%==============================|
% Plots of the Relative Errors |
%==============================|
figure(1); 
shift = 1/(H(end)*2);
subplot(1,2,1);
loglog(H,S_err,'o-b','LineWidth',2);         hold on 
loglog(H,H*S_err(end)*shift,':b','LineWidth',1.5); 
xlabel('h');    ylabel('E(h)');    axis square tight
title('Relative Error on S(t)');
xticks(10.^(-5:-1)); yticks(10.^(-10:-1));
subplot(1,2,2); 
loglog(H,phi_err,'s-c','LineWidth',2); hold on;
loglog(H,H*shift*phi_err(end),':c','LineWidth',1.5); 
xlabel('h');     ylabel('E(h)');   axis square tight
title('Relative Error on φ(t)');
xticks(10.^(-5:-1)); yticks(10.^(-10:-1));

fprintf('\nDone.\n\n');

% ==============================================================================
% End of Test1_Paper_JCD.m
% ==============================================================================
