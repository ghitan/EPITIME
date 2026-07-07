% Test 2: Unimodal memory kernel (part of EPITIME software)
%
% B. Buonomo, E. Messina, C. Panico, A. Vecchio.
% "An integral renewal equation approach to behavioural epidemic models with 
%  information index"
% Journal of Mathematical Biology, 2025, 90(8)
% DOI : 10.1007/s00285-024-02172-y
%
% ==============================================================================
% Version:    1.0
% Date:       April, 2026
% Authors:    C. Panico, G. Zanghirati
% See:        https://arxiv.org/abs/2605.00067v1
% License:    GNU GPL v.3
% Repository: https://github.com/ghitan/EPITIME

clear; close all; clc
fprintf('\nEPITIME behavioural: test2: unimodal memory kernel\n');

%=========================================================
%% Input Parameters 
T  = 3000;
h  = 0.1;
N  = 5*(10^7);
mu = 1.0/(75*365); 

%% Functions:
    %% Infectivity function: 
    % In order to consider the trapezoidal infectivity function, we follow the approach proposed in [1]. 
    % We assume that the basic reproduction number is equal to 3.3, representing SARS, 
    % and this value is available in the literature. If R0 = 3.3, then the function A is given as follows, 
    % where the parameter p0 allows us to satisfy the definition of R0 through the integral of A.
    R0 = 3.3; 
    R0 = R0/N;
    tp = linspace(4,14);
    for i = 1:100
        A(i) = kernel_sars(tp(i));
    end
    Ahat0 = trapz(tp,A);
    p0 = R0/(Ahat0);
    A  = @(x)( p0.*kernel_sars(x) );

    %% Memory kernel
    a = 1/30;
    K = @(x)( (a^2).* x .*exp(-a.*x) );

    %% Message function
    g = @(x)( x );

    %% Inhibition function
    alpha = 5*(10^4);
    betaM = @(x)( 1./(1+alpha*x) ); % Represents the inhibition function used in [1] and [2] 

%% Equilibria given in the standardized form
R0n = R0*N;
Feq = mu/N*(R0n-1)/(1+alpha*mu); 
lam = mu;
Seq = lam*(1/N+alpha*Feq)/(mu*(R0n/N));

%% Initial data
S0  = 0.7*N;

fprintf('\nCalling the solver...');
%=========================================================   
% Format: NSFD_behavioural(T, h, N, mu, S0, A, K, betaM, g, compGamh, Verbosity)
[t, Y, P] = NSFD_behavioural(T, h, N, mu, S0, A, K, betaM, g, 1, 1);
%=========================================================
fprintf('done.\n');

%=========================================================
%                 ----- FIGURES -----
% Variable extraction
s = Y(1,:); 
f = Y(2,:);
  
plotData(t, f, Feq)

fprintf('\nDone.\n\n');

%=========================================================
%       ---- TRAPEZOIDAL INFECTIVITY FUNCTITON ------
%                                                                                                                  
function A = kernel_sars(tau)

    ta = 4; tb = 7; tc = 11; td = 14;  
    mu = 1 / (75 * 365);

    A = zeros(size(tau));

    mask1 = (tau > ta) & (tau < tb);
    mask2 = (tau >= tb) & (tau <= tc);
    mask3 = (tau > tc) & (tau < td);

    A(mask1) = exp(-mu .* tau(mask1)) .* (tau(mask1) - ta) ./ (tb - ta);
    A(mask2) = exp(-mu .* tau(mask2));
    A(mask3) = exp(-mu .* tau(mask3)) .* (td - tau(mask3)) ./ (td - tc);
end

% ==============================================================================
% End of Behavioural_Trapezoidal_IF.m
% ==============================================================================