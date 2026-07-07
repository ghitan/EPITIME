% TEST 1 of the Paper: AGE OF INFECTION WITH LOW REGULARITY KERNEL (part of EPITIME software)
%
% E. Messina, M. Pezzella, A. Vecchio.
% "Nonlocal finite difference discretization of a class of renewal equation
%  models for epidemics"
% Mathematical Biosciences and Engineering, 2023, 20(7): 11656-11675.
% DOI: 10.3934/mbe.2023518
%
% For further details on the Age of Infection model we refer also to:
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

close all; clear; clc
fprintf('\nEPITIME general: test 1: AoI with a low-regularity kernel\n');

%====================|
% Problem definition |
%====================|
% Two alternative input modes
% (1) Structured input:
% => NSFD_Renewal(Ren_Problem)
N = 5e+4;
Ren_Problem       = struct();
Ren_Problem.S0    = 49950;
Ren_Problem.A     = @(i,t)(  ((t > 14) & (t < 16)) .* ((t - 14).*(19 - t)/14) ...
                           + ((t >= 16) & (t < 19)) .* ((19 - t)/7) );
Ren_Problem.phi0  = @(i,t)( (N - Ren_Problem.S0) * Ren_Problem.A(i,t) );
Ren_Problem.P0    = @(t)( zeros(size(t)) );
Ren_Problem.B     = @(t)( zeros(size(t)) );
Ren_Problem.T     = 1e+2;
Ren_Problem.h     = 1e-2;
Ren_Problem.beta  = 1e-4;
Ren_Problem.BETA  = 1;
Ren_Problem.alpha = 0;
Ren_Problem.c     = 0;
Ren_Problem.verbosity = 0;

%%% (2) Positional input:
%%% => NSFD_Renewal(S0, phi0, P0, A, B, beta, BETA, alpha, c, h, T)
% A = @(i,t) ((t > 14) & (t < 16)) .* ((t - 14).*(19 - t)/14) ...
%    +((t >= 16) & (t < 19)) .* ((19 - t)/7);
% beta = 1e-4;       N = 5e+4;     S0 = 49950;      h = 1e-2;
% T = 1e+2;       phi0 = @(n,t) (N-S0)*A(n,t);      P0 = @(t) zeros(size(t)); 
% B = @(t) zeros(size(t));        BETA = 1;         alpha = 0;     c = 0;

fprintf('\nComputing the solution...')
[t, S, phi, P] = NSFD_Renewal( Ren_Problem );
fprintf('done.\n')

%==============================|
% Continuous R₀ and Final Size |
%==============================|
tmpA = @(t)( Ren_Problem.A(1,t) );
R0 = Ren_Problem.beta * N * integral(tmpA,0,inf);
FS_Function = @(x)( log(Ren_Problem.S0/x) - R0 * (1 - x/N) );
options = optimset('TolX',1e-14,'TolFun',1e-14);
S_infty = fzero(FS_Function, [1e-14, Ren_Problem.S0], options);
S_inf_num = S(end);
err_infty = abs(S_inf_num - S_infty) / S_infty;
fprintf('\nThe Continuous Basic Reproduction Number is R0 = %.2f \n', R0); 
fprintf('The Continuous Final Size is S(∞) = %.5e \n', S_infty); 
fprintf(['The Numerical Final Size is S_∞ = %.5e with a relative error of ' ...
         '%.5e \n'], S_inf_num, err_infty); 

%================================|
% Plots of the Numerical Soution |
%================================|
figure(1);
plot(t, S, 'b', 'Linewidth',2); hold on; plot(t, phi, 'c', 'Linewidth',2);
scatter(t(end), S_infty, 'LineWidth', 2, ...
        MarkerEdgeColor = 'b',MarkerFaceColor = 'b');
legend('S(t)','φ(t)','S(∞)'); xlabel('t')
title(sprintf('NSFD Numerical Solution with h = %.2e', Ren_Problem.h))

fprintf('\nDone.\n\n');

% ==============================================================================
% End of Test1_Paper_MBE.m
% ==============================================================================