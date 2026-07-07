% BEHAVIOURAL EPIDEMIC PEAKS INVESTIGATION (part of EPITIME software)
%
% ==========================================================================
% Version:    1.0
% Date:       June, 2026
% Authors:    C. Panico, G. Zanghirati
% See:        https://arxiv.org/abs/2605.00067v1
% License:    GNU GPL v.3
% Repository: https://github.com/ghitan/EPITIME

close all; clear; clc
fprintf('\nEPITIME behavioural: epidemic peaks investigation\n');

%=========================================================
%% Input Parameters 
T  = 1000;
h  = 0.1;
N  = 5*(10^7);
mu = 1.0/(75*365); 

%% Functions:
ni = 1/7; 
R0 = 20; 
bb = (R0/N)*(mu+ni)^2; 
A  = @(x)( bb.*x.*exp(-(mu+ni).*x) );

a = 1/30;
K = @(x) (a).*exp(-a.*x);

g = @(x) x;

alpha = 8*(10^3);
betaM = @(x) 1./(1+alpha*x); 

%% Equilibria 
Feq = mu/N*(R0-1)/(1+alpha*mu); 
lam = mu;
Seq = lam*(1/N+alpha*Feq)/(mu*(R0/N));

%% Initial data
S0  = 0.99*Seq*N;

fprintf('\nPeaks estimation started...');
%=========================================================   
% Solver execution
[t, Y, P] = NSFD_behavioural(T, h, N, mu, S0, A, K, betaM, g, 1, 1);
%=========================================================

%% Variable Extraction 
S_N = Y(1,:) * N; 
F_N = Y(2,:) * N; 

Nt = length(t);
incidence_NS = zeros(1, Nt-1);
for n = 1:Nt-1
    D = h;
    if (T - floor(T/h)*h ~= 0 && n == Nt-1)
        D = T - floor(T/h)*h;
    end
    incidence_NS(n) = N*((Y(1,n)-Y(1,n+1))/h+mu*(1-Y(1,n+1)));
end
incidence_NS(Nt) = incidence_NS(Nt-1); 

%=========================================================
%%  NUMERICAL PEAK DETECTION 
%=========================================================

% 1. Peaks of the Driving Factor F(t)
f_diff1 = F_N(2:end-1) > F_N(1:end-2);
f_diff2 = F_N(2:end-1) > F_N(3:end);
idx_f = find(f_diff1 & f_diff2) + 1;

pks_f = F_N(idx_f);
locs_f = t(idx_f);
valid_f = pks_f > (max(F_N) * 0.01);
pks_f = pks_f(valid_f);
locs_f = locs_f(valid_f);

% 2. Peaks of the Incidence
inc_diff1 = incidence_NS(2:end-1) > incidence_NS(1:end-2);
inc_diff2 = incidence_NS(2:end-1) > incidence_NS(3:end);
idx_inc = find(inc_diff1 & inc_diff2) + 1;

pks_inc = incidence_NS(idx_inc);
locs_inc = t(idx_inc);
valid_inc = pks_inc > (max(incidence_NS) * 0.01);
pks_inc = pks_inc(valid_inc);
locs_inc = locs_inc(valid_inc);

fprintf('\nPeaks estimation finished.');

%% Output
fprintf('\n=========================================================\n');
fprintf('  RESULTS \n');
fprintf('=========================================================\n');
fprintf('\n--- PEAKS OF DRIVING FACTOR F(t) ---\n');
if isempty(locs_f)
    fprintf('No peaks detected.\n');
else
    for i = 1:length(locs_f)
        fprintf('Peak %d: Time = %6.2f, F(t) = %12.4e\n', i, locs_f(i), pks_f(i));
    end
end

fprintf('\n--- PEAKS OF INCIDENCE ---\n');
if isempty(locs_inc)
    fprintf('No peaks detected.\n');
else
    for i = 1:length(locs_inc)
        fprintf('Peak %d: Time = %6.2f, Incidence(t) = %12.2f\n', i, locs_inc(i), pks_inc(i));
    end
end
fprintf('=========================================================\n');

%% Data Visualization
figure('Name', 'Epidemic Peaks', 'NumberTitle', 'off');

% Subplot 1: Driving Factor F(t)
subplot(2,1,1);
plot(t, F_N, 'b-', 'LineWidth', 1.5); hold on;
plot(locs_f, pks_f, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
grid on;
%title('Peak Detection in Driving Factor F(t) ');
xlabel('t');
ylabel('F(t)');
legend('F(t)', 'Detected Peaks', 'Location', 'best');

% Subplot 2: Incidence
subplot(2,1,2);
plot(t, incidence_NS, 'g-', 'LineWidth', 1.5); hold on;
plot(locs_inc, pks_inc, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
grid on;
%title('Peak Detection in Incidence');
xlabel('t');
ylabel('Inc(t)');
legend('Inc(t)', 'Detected Peaks', 'Location', 'best');

fh = gcf;

% -----------------------------
% Vector graphics export
% -----------------------------
fname = 'peaks_beh.pdf';
exportgraphics(fh, fname, ...
               'BackgroundColor', 'none', ...
               'Resolution', 1200, ...
               'ContentType', 'vector');

fprintf('\n\nDone.\n\n');

% ==============================================================================
% End of Behavioural_Peaks_Estmation.m
% ==============================================================================