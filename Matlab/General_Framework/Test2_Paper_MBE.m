% TEST 2 of the Paper: AOI WITH SYMPTOMATIC AND ASYMPTOMATIC INFECTIONS (part of EPITIME software)
%
% E. Messina, M. Pezzella, A. Vecchio.
% "Nonlocal finite difference discretization of a class of renewal equation
%  models for epidemics"
% Mathematical Biosciences and Engineering, 2023, 20(7): 11656-11675.
% DOI: 10.3934/mbe.2023518
%
% For further details on the Age of Infection model with with symptomatic 
% and asymptomatic infections we refer also to:
% Bai F.
% "An age-of-infection model with both symptomatic and asymptomatic
%  infections"
% Journal of Mathematical Biology, 2023, 86(5): 239–252.
% DOI: 10.1007/s00285-023-01920-w
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
% Two alternative input modes
% (1) Structured input:
% => NSFD_Renewal( Ren_Problem )
N = 1e+4;         a = 1;            f = 0.783;
Ren_Problem       = struct();
Ren_Problem.S0    = [0.9*N; 0.9*N]; 
Ren_Problem.A     = @(i,t)(  (i == 1) .* (5*gampdf(t, 1, 2)) .* exp(-0.5*sqrt(t)) ...
                           + (i == 2) .* (1./(1 + 0.6*t))  .* (gampdf(t, 3, 2)) );
Ren_Problem.phi0  = @(i,t)( (N - Ren_Problem.S0(1)) * Ren_Problem.A(i,t) / ...
                            ( f*(i == 1) + (1-f)*(i == 2) ) ); 
Ren_Problem.P0    = @(t)( zeros(size(t)) );
Ren_Problem.B     = @(t)( zeros(size(t)) ); 
Ren_Problem.T     = 30;
Ren_Problem.h     = 1e-2;
Ren_Problem.beta  = [a/N; a/N];
Ren_Problem.BETA  = ones(2,2);
Ren_Problem.alpha = zeros(2,1);
Ren_Problem.c     = zeros(2,1);
Ren_Problem.verbosity = 0;

%%% (2) Positional input:
%%% => NSFD_Renewal(S0, phi0, P0, A, B, beta, BETA, alpha, c, h, T)
% N = 1e+4;     S0=[0.9*N; 0.9*N];      a=1;
% beta = [a/N; a/N];        h=1e-2;     T=30;       f=0.783;
% A = @(i,t)(  (i == 1) .* (5*gampdf(t, 1, 2)) .* exp(-0.5*sqrt(t)) ...
%            + (i == 2) .*(1./(1 + 0.6*t))  .* gampdf(t, 3, 2)  );
% phi0 = @(n,t)( (N - S0(1)) * A(n,t) / (f*(n==1) + (1-f)*(n==2)) );
% P0 = @(t)( zeros(size(t)) ); B = @(t)( zeros(size(t)) );
% BETA = ones(2,2);    alpha = zeros(2,1);     c = zeros(2,1);

[t, S, phi, P, Perf] = NSFD_Renewal( Ren_Problem );

%=============================================================|
% Number of symptomatic and asymptomatic infected individuals |
%=============================================================|
Bas  = @(i,t)( (i==1) .* exp(-0.5*sqrt(t)) + (i==2) .* (1./(1 + 0.6*t)) );
I0   = @(t)( f*(N-Ren_Problem.S0(1))*Bas(1,t) );           I0vals = I0(t);
A0   = @(t)( (1-f)*(N-Ren_Problem.S0(1))*Bas(2,t) );       A0vals = A0(t);
Nt   = size(S,2);
pi1  = 0.5;                                  pi2 = pi1;
I    = zeros(Nt,1);                           AI = zeros(Nt,1);
I(1) = I0(0);                              AI(1) = A0(0);
Bas_vals = zeros(2,Nt);
phi_sum = sum(phi,1);
for i = 1:2
    Bas_vals(i,:) = Bas(i,t);
end
tmph = Ren_Problem.h * a/N;
for n = 1:Nt-1
    I(n+1)  = I0vals(n) + tmph*f * sum(S(1,2:n+1) .* phi_sum(1:n) ...
              .* Bas_vals(1,n:-1:1), 2);
    AI(n+1) = A0vals(n) + tmph*(1-f) * sum(S(2,2:n+1) .* phi_sum(1:n) ...
              .* Bas_vals(2,n:-1:1), 2);
end

%=====================================|
% Continuous and Numerical Final Size |
%=====================================|
F = a*integral( @(t) Ren_Problem.A(1,t), 0, inf ) ...
    + a*integral(@(t) Ren_Problem.A(2,t) , 0, inf);
s0 = Ren_Problem.S0(1);
FS_Function = @(x) log(s0/x) - F*( (s0-x)/N) - a/N*(integral(@(t) ...
    Ren_Problem.phi0(1,t),0,inf) + integral( @(t) Ren_Problem.phi0(2,t), ...
    0, inf) );
options = optimset('TolX',1e-14,'TolFun',1e-14);
S_infty = fzero(FS_Function, [1e-14, s0], options);
S_inf_num = S(end);
err_infty = abs(S_inf_num - S_infty) / S_infty;
fprintf('The Value of F is %.2f \n', F); 
fprintf('The Continuous Final Size is S(∞) = %.5e \n', S_infty); 
fprintf(['The Numerical Final Size is S_∞ = %.5e with a relative error of ' ...
         '%.5e \n'], S_inf_num, err_infty); 

%================================|
% Plots of the Numerical Soution |
%================================|
figure()
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, S(1,:), 'b', 'LineWidth', 2); hold on;
plot(t, N - (S(1,:)' + I + AI), 'c', 'LineWidth', 2);
scatter(t(end), S_infty, 'filled', ...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
legend('S(t)', 'R(t)', 'S(∞)', 'Location', 'southeast');
xlabel('t');            axis square;
nexttile;
plot(t, I, 'Color', 'magenta', 'LineWidth', 2); hold on;
plot(t, AI, 'r', 'LineWidth', 2);
legend('Iˢ(t)', 'Iᵃ(t)');
xlabel('t');            axis tight square;
sgt = sgtitle(sprintf('NSFD Numerical Solution with h = %.2e', Ren_Problem.h));
sgt.FontSize = 14;
sgt.FontWeight = 'bold';

