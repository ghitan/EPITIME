% Simple example of an Age of Infection model (part of EPITIME software)
%
%========================================================================%
%                          PROBLEM DEFINITION                            %
%------------------------------------------------------------------------%
% TWO ALTERNATIVE INPUT MODES:                                           %
% 1. OPTION 1: Structured input (using the AoI_problem struct)           %
% 2. OPTION 2: Positional input (using individual variables)             %
%                                                                        %
% INSTRUCTIONS: Please comment out the option that you are NOT using.    %
%                                                                        %
% TIP: To comment/uncomment multiple lines quickly:                      %
%      - Windows: Select lines + [Ctrl + R] to Comment, [Ctrl + T] to Un %
%      - Mac:     Select lines + [Cmd  + R] to Comment, [Cmd  + T] to Un %
%========================================================================%
% 
% Version:    1.0
% Date:       April, 2026
% Authors:    M. Pezzella, G. Zanghirati
% See:        https://arxiv.org/abs/2605.00067v1
% License:    GNU GPL v.3
% Repository: https://github.com/ghitan/EPITIME

clearvars; clc; close all;

% OPTION 1: STRUCTURED INPUT_________________________________________________
% AoI_problem = struct();                                                  %|
% AoI_problem.T = 60;                                                      %|
% AoI_problem.h = 5e-3;                                                    %|
% AoI_problem.N = 1e3;                                                     %|
% AoI_problem.S0 = 0.99 * AoI_problem.N;                                   %|
% lambda = 1; k = 9;                                                       %|
% AoI_problem.A = @(t)( 0.3*t.^(k-1) .* exp( -t/lambda ) / ( lambda.^k ... %|
%     *factorial(k-1) ) );                     % Erlang distribution kernel |
% AoI_problem.phi0 = @(t)( (AoI_problem.N - AoI_problem.S0) ...            %|
%     * AoI_problem.A(t) );                                                %|
% AoI_problem.beta = 5e-3;                                                 %|
% AoI_problem.verbosity = 1;                                               %|
% [t, Y, P] = NSFD_AoI( AoI_problem );                                     %|
% __________________________________________________________________________|

% OPTION 2: POSITIONAL INPUT:______________________________________________
T = 60;                                                                  %|
h = 5e-3;                                                                %|
N = 1e3;                                                                 %|
S0 = 0.99 * N;                                                           %|
lambda = 1; k = 9;                                                       %|   
A = @(t)( 0.3*t.^(k-1) .* exp( -t/lambda ) / ( lambda.^k*factorial(k-1)) );
phi0 = @(t)( (N - S0) * A(t) );                                          %|
beta = 5e-3;                                                             %|
verbosity = 1;                                                           %|
[t, Y, P]= NSFD_AoI(S0, N, T, h, beta, A, phi0, verbosity);              %|
%_________________________________________________________________________|

%===============================================|
% Plots of the solutions (For both the options) |
%===============================================|
figure()
plot(t,Y(1,:),"LineWidth",2);
xlabel('t'); ylabel('S(t)')
title('Susceptible Individuals')

figure()
plot(t,Y(2,:),"LineWidth",2);
xlabel('t'); ylabel('φ(t)')
title('Mean Infectivity')

%===============================================|
% Plots of the kernel function                  |
%===============================================|
if exist('AoI_problem', 'var') && isfield(AoI_problem, 'A')
    figure()
    plot(t, AoI_problem.A(t), "LineWidth", 2);
    xlabel('t'); ylabel('A(t)')
    title('Infectivity Kernel')
elseif exist('A', 'var') || exist('A', 'file')
    figure()
    plot(t, A(t), "LineWidth", 2); 
    xlabel('t'); ylabel('A(t)')
    title('Infectivity Kernel')
else
    disp('Since the default kernel is employed, no plotting for A is shown.');
end
