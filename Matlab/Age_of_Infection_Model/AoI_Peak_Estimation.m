% EPIDEMIC PEAK INVESTIGATION (part of EPITIME software)
%
% ==========================================================================
% Version:    1.0
% Date:       June, 2026
% Authors:    M. Pezzella, G. Zanghirati
% See:        https://arxiv.org/abs/2605.00067v1
% License:    GNU GPL v.3
% Repository: https://github.com/ghitan/EPITIME

clear; close all; clc
fprintf('\nEPITIME AoI: epidemic peak investigation\n');

N = 100 ;
i0 = [0.0268, 1.27e-6, 0.05, 1.33e-4, 5e-6, 0.001];
s0 = 1-i0;
betaN = [4.6291, 0.5, 10, 0.13905, 0.3333, 0.5]; 
gamma = [2.82, 0.3, 1, 0.018379, 0.1111, 0.3];
N_cases = length(i0);
peak_Continuous = zeros(N_cases,1);
t_peak = zeros(N_cases,1);
t_peak_Sir = zeros(N_cases,1);
S_peak = zeros(N_cases,1); E_peak = S_peak;

% Prepare the running echo to console
Ndisp = max(1., N_cases / 10.) * ( 1 : min(N_cases, 10) )';
k = 1;
fprintf('\nPeak estimation started...\n0%%..');

for j = 1:N_cases
    AoI_problem.S0 = s0(j)*N;
    AoI_problem.N = N;
    AoI_problem.T = 100;
    AoI_problem.h = 1e-3;
    AoI_problem.beta = betaN(j)/N;
    AoI_problem.A = @(t) exp(- gamma(j) * t);
    AoI_problem.phi0 = @(t) AoI_problem.N * i0(j) * AoI_problem.A(t);
    AoI_problem.verbosity = 1;
    R0 = AoI_problem.beta * AoI_problem.N * integral(AoI_problem.A, 0, inf);
    
    s_tmp = s0(j);
    i_tmp = i0(j);
    sm = 1/R0;
    peak_integr = @(s) 1/betaN(j) .* 1./( s .* (s-sm .* log(s./s_tmp) ...
                       -i_tmp - s_tmp ));
    peak_Continuous(j) = integral(peak_integr, s_tmp, sm, ...
                                  'RelTol',1e-12,'AbsTol',1e-12);
    
    [t, Y, P] = NSFD_AoI(AoI_problem);
    incidence = AoI_problem.beta * Y(1,:) .* Y(2,:);
    
    Sir = @(t,y) [-betaN(j)*y(1)*y(2); betaN(j)*y(1)*y(2)-gamma(j)*y(2);...
                  gamma(j)*y(2)];
    opts = odeset('RelTol',1e-13,'AbsTol',1e-13);
    [t_Sir, Y_Sir] = ode89(Sir, [0, AoI_problem.T], [s_tmp; i_tmp; 1-s_tmp-i_tmp], opts);
    
    [~, idx_inc] = max(incidence);
    [~, idx] = min(abs(Y(1,:)/N - 1/R0));
    t_peak(j) = t(idx);
    S_peak(j) = Y(1,idx);
    E_peak(j) = AoI_problem.beta * Y(2,idx);
    [~, idx_Sir] = min(abs(Y_Sir(:,1) - 1/R0));
    t_peak_Sir(j) = t_Sir(idx_Sir);
    % Display running iterations percentage
    if ( j >= Ndisp(k) ), fprintf('%d%%..', fix(100. * Ndisp(k)/ N_cases)); k = k + 1; end
end
fprintf('%c%c\nPeak estimation finished.\n', char(8), char(8));

%
format short e

fprintf('\n');
fprintf('-----------------------------------------------------------------------------------------------------------------\n');
fprintf('    s0         i0       beta     gamma         t_p*     t_p(AoI_NSFD)  R.E.(t_p)       E_p\n');
fprintf('-----------------------------------------------------------------------------------------------------------------\n');

for j = 1:N_cases
    fprintf('%9.5f  %9.2e  %8.5f  %8.5f  %11.5f  %11.5f  %11.2e  %11.5f\n', ...
        s0(j), ...
        i0(j), ...
        betaN(j), ...
        gamma(j), ...
        peak_Continuous(j), ...
        t_peak(j), ...
        abs(t_peak(j)-peak_Continuous(j))/peak_Continuous(j), ...
        E_peak(j));
end

fprintf('-----------------------------------------------------------------------------------------------------------------\n');


%% ================================================================
% Convergence study for the epidemic peak time
% ================================================================
h_vec = [1e-1, 5e-2, 2.5e-2, 1.25e-2, 6.25e-3, 3.125e-3, 1.5625e-3];
Nh = length(h_vec);
t_peak_h = zeros(Nh,N_cases);
err_tp   = zeros(Nh,N_cases);
eoc_tp   = NaN(Nh,N_cases);

% Prepare the running echo to console
% Ndisp = max(1., N_cases / 10.) * ( 1 : min(N_cases, 10) )'; 
kk = 1;
fprintf('\nEPITIME AoI: convergence study started...\n0%%..');

for j = 1:N_cases
    % Problem definition
    AoI_problem.S0 = s0(j)*N;
    AoI_problem.N = N;
    AoI_problem.T = 100;
    AoI_problem.beta = betaN(j)/N;
    AoI_problem.A = @(t) exp(-gamma(j)*t);
    AoI_problem.phi0 = @(t) N*i0(j)*AoI_problem.A(t);
    AoI_problem.verbosity = 0;
    R0 = AoI_problem.beta*AoI_problem.N*integral(AoI_problem.A,0,inf);
    % Exact peak time
    tp_exact = peak_Continuous(j);
    for k = 1:Nh
        AoI_problem.h = h_vec(k);
        [t,Y,P] = NSFD_AoI(AoI_problem);
        [~,idx] = min(abs(Y(1,:)/N - 1/R0));
        t_peak_h(k,j) = t(idx);
        err_tp(k,j) = abs(t_peak_h(k,j)-tp_exact)/tp_exact;
        if k > 1
            eoc_tp(k,j) = log(err_tp(k-1,j)/err_tp(k,j))/log(2);
        end
    end
    % Display running iterations percentage
    if ( j >= Ndisp(kk) ), fprintf('%d%%..', fix(100. * Ndisp(kk)/ N_cases)); kk = kk + 1; end
end
fprintf('%c%c\nConvergence study finished.\n', char(8), char(8));

format short e
for j = 1:N_cases
    fprintf('\n');
    fprintf('===============================================================\n');
    fprintf('Case %d: i0 = %.2e, beta = %.5f, gamma = %.5f\n',...
            j,i0(j),betaN(j),gamma(j));
    fprintf('===============================================================\n');
    fprintf('      h           t_p(h)         Rel. Error         EOC\n');
    fprintf('---------------------------------------------------------------\n');
    for k = 1:Nh
        if k == 1
            fprintf('%10.2e   %12.6f    %12.2e        ----\n',...
                    h_vec(k),t_peak_h(k,j),err_tp(k,j));
        else
            fprintf('%10.2e   %12.6f    %12.2e      %6.4f\n',...
                    h_vec(k),t_peak_h(k,j),err_tp(k,j),eoc_tp(k,j));
        end
    end
    fprintf('===============================================================\n');
end


%% =========================================================================
%  GRAPHICS SETTINGS
%
%  Define font sizes and font family used throughout all figures.
% =========================================================================
axisfontsize   = 11;
labelfontsize  = 11;
titlefontsize  = 13;
legendfontsize = labelfontsize;

arcstr = computer('arch'); OS = lower( arcstr(1) );
switch OS
    case 'w', fontname = 'Palatino Linotype';   % Windows
    case 'm', fontname = 'Palatino';            % macOS
    case 'g', fontname = 'URW Palladio L';      % Linux (GLNX)
    otherwise, warning('Unknown underlying operating system: %s', arcstr);
end

% ================================================================
%  Convergence plot: relative error of peak time
% ================================================================
figure;
colors = jet(N_cases);
for j = 1:N_cases
    loglog(h_vec, err_tp(:,j), '-o', ...
           'LineWidth', 1.5, ...
           'MarkerSize', 5, ...
           'Color', colors(j,:)); hold on;
end
% Reference slope O(h)
h_ref = h_vec;
C = 0.5* err_tp(end,1) / h_vec(end);   % just a scaling for visual guide
loglog(h_ref, C*h_ref, 'k--', 'LineWidth', 1.5);
xlabel('h', 'FontName', fontname, 'FontSize', labelfontsize);
ylabel('Relative error on $t_p$', ...
       'Interpreter','latex', ...
       'FontName', fontname, ...
       'FontSize', labelfontsize);
legend_entries = cell(1,N_cases);
for j = 1:N_cases
    legend_entries{j} = sprintf('Case %d', j);
end
legend([legend_entries, {'$\mathcal{O}(h)$'}], ...
       'Interpreter','latex', ...
       'FontSize', legendfontsize, ...
       'Location','northwest');
set(gca, ...
    'FontName', fontname, ...
    'FontSize', axisfontsize);
title('AoI\_NSFD - Convergence of epidemic peak time', ...
      'FontName', fontname, ...
      'FontSize', titlefontsize);
hold off;

axh = gca;
set(axh,'FontName',fontname,'FontSize',axisfontsize);
width_cm = 16.51;
width_in = width_cm / 2.54;
fh = gcf;
set(fh, 'Units','inches');
pos = get(fh, 'Position');        % [left bottom width height]
aspect = pos(4) / pos(3);         % height / width
new_width  = width_in;
new_height = width_in * aspect;
set(fh, 'Position', [pos(1) pos(2) new_width new_height]);

fprintf('\n\nDone.\n\n');

% ==============================================================================
% End of AoI_Peak_Estmation.m
% ==============================================================================