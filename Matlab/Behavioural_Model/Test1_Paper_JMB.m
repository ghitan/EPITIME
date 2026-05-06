% Test 1: Exponential memory kernel (part of EPITIME software)
% 
% B. Buonomo, E. Messina, C. Panico, A. Vecchio.
% "An integral renewal equation approach to behavioural epidemic models with 
%  information index"
% Journal of Mathematical Biology, 2025, 90(8).
% DOI : 10.1007/s00285-024-02172-y
%
% ==========================================================================
% Version:    1.0
% Date:       April, 2026
% Authors:    C. Panico, G. Zanghirati
% See:        https://arxiv.org/abs/2605.00067v1
% License:    GNU GPL v.3
% Repository: https://github.com/ghitan/EPITIME

clear; clc

%=========================================================
%% Input Parameters 
T  = 1000;
h  = 0.1;
N  = 5*(10^7);
mu = 1.0/(75*365); 

%% Functions:
    %% Infectivity function: 
    % In order to consider the unimodal infectivity function, we follow the approach proposed in [1]. 
    % We assume that the basic reproduction number is equal to 20, representing measles, 
    % and this value is available in the literature. If R0 = 20, then the function A is given as follows, 
    % where the parameter bb allows us to satisfy the definition of R0 through the integral of A.
    ni = 1/7; 
    R0 = 20; 
    bb = (R0/N)*(mu+ni)^2; 
    A  = @(x)( bb.*x.*exp(-(mu+ni).*x) );

    %% Memory kernel
    a = 1/30;
    K = @(x) (a).*exp(-a.*x);

    %% Message function
    g = @(x) x;

    %% Inhibition function
    alpha = 8*(10^3);
    betaM = @(x) 1./(1+alpha*x); % Represents the inhibition function used in [1] 

%% Equilibria given in the standardized form
Feq = mu/N*(R0-1)/(1+alpha*mu); 
lam = mu;
Seq = lam*(1/N+alpha*Feq)/(mu*(R0/N));

%% Initial data
S0  = 0.99*Seq*N;

%=========================================================   
% Format: NSFD_behavioural(T, h, N, mu, S0, A, K, betaM, g, Verbosity)
[t, Y, P] = NSFD_behavioural(T, h, N, mu, S0, A, [], betaM, g, 1);
%=========================================================


%=========================================================
%                    -----FIGURES-----
% Variable extraction
s = Y(1,:); 
f = Y(2,:);

plotData(t, f, Feq)
