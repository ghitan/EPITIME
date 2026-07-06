% TEST 2 of the Paper: LONG TIME BEHAVIOUR INVESTIGATION (part of EPITIME software)
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

close all; clear; clc
fprintf('\nEPITIME AoI: test 2: long time behaviour investigation\n');

%====================|
% Problem definition |
%====================|
% Two alternative input modes
% (1) Structured input:
% => NSFD_AoI( AoI_problem )
mu = 0.2;              sigma = 2*mu;
AoI_problem = struct();
AoI_problem.T = 30;
AoI_problem.h = 1e-1;
AoI_problem.A = @(t) ( 1 / (sigma * sqrt(2 * pi)) ) *...
                     exp( -0.5 * ((t - mu) ./ sigma).^2 ); % ORIGINAL PAPER TEST
                                     % OTHER TESTS: 
                                     % AoI_problem.A = @(t)( exp(-t) );
                                     % AoI_problem.A = @(t)( 1./(1+t.^4) );
AoI_problem.N = 1e5;
AoI_problem.beta = 3e-5;
AoI_problem.S0 = 99950;
AoI_problem.phi0 = @(t)( (AoI_problem.N - AoI_problem.S0) * AoI_problem.A(t) );
AoI_problem.verbosity = 0;
% (2) Positional input:
%%% => NSFD_AoI(S0, N, T, h, beta, A, phi0, verbosity)
% T = 30;  N = 1e5;   S0 = 99950;   beta = 3e-5;   h=1e-1;  verbosity = 0;
% A = @(t)( normpdf(t,mu,sigma) ); phi0 = @(t)( (N - S0) .* A(t) );

TestN = 3;                   
S_final = zeros(TestN,1);                 phi_final = zeros(TestN,1);  
H = zeros(TestN,1);                       R0_num = zeros(TestN,1);
betaN = AoI_problem.beta * AoI_problem.N;
for j = 1:TestN
    h_test = 10^(-j);
    H(j) = h_test;
    AoI_problem_test = AoI_problem;            AoI_problem_test.h = h_test;
    [~, Y_test, ~] = NSFD_AoI( AoI_problem_test );
    % [~, Y_test, ~] = NSFD_AoI(S0, N, T, h_test, beta, A, phi0);  % Alt. (2)
    S_final(j) = Y_test(1,end);                phi_final(j) = Y_test(2,end);
    % S_final = S_∞(h) is the Numerical Final Size
    R0_num(j) = betaN * h_test * sum( AoI_problem.A(h_test:h_test:1e+4) );
end

%==============================|
% Continuous R₀ and Final Size |
%==============================|
R0 = betaN * integral(AoI_problem.A, 0, inf);
FS_Function = @(x)( log(AoI_problem.S0/x) - R0 * (1 - x/AoI_problem.N) );
options = optimset('TolX',1e-14,'TolFun',1e-14);
S_infty = fzero(FS_Function, [1e-14, AoI_problem.S0], options);
err_infty = abs(S_final - S_infty) / S_infty;
fprintf('\nThe Continuous Basic Reproduction Number is R0 = %.2f \n', R0); 
fprintf('The Continuous Final Size is S(∞) = %.5e \n', S_infty); 

%================================================================|
% Table with Num. R₀, Num. Final Size and Num. Final Infectivity |
%================================================================|
fprintf('\n %-9s %-15s %-14s %-17s %-15s\n', ...
        'h', 'R₀(h)', 'S_∞(h)', 'Rel. S_∞ Err. ', 'φ_∞(h)');
for j = 1:length(H)
    fprintf('10^-%-5d %-15.5e %-15.5e %-18.5e %-15.5e\n', ...
            j, R0_num(j), S_final(j), err_infty(j), phi_final(j));
end

% =========================================================================
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
width_cm = 16.51;
width_in = width_cm / 2.54;

%============================================|
% Plot of the numerical solution  for h=1e-3 |
%============================================|
AoI_problem_plot = AoI_problem_test;             AoI_problem_plot.h = 1e-3;
[t_num, Sol_num, ~] = NSFD_AoI( AoI_problem_plot );
figure()
plot(t_num,Sol_num(1,:),'LineWidth',2,'Color','b'); hold on
plot(t_num,Sol_num(2,:),'LineWidth',2, 'Color','c');
yline(S_infty,'--', 'Color',[0.00,0.74,0.78],'LineWidth', 2);
xlh = xlabel('t'); ylh = ylabel('');
lh = legend('S(t)','φ(t)','S(∞)');
th = title(sprintf('NSFD Numerical Solution with h = %.2e', AoI_problem_plot.h));
axh = gca;
set(axh,'FontName',fontname,'FontSize',axisfontsize);
set([xlh,ylh],'FontName',fontname,'FontSize',labelfontsize);
set(th,'FontName',fontname,'FontSize',titlefontsize,...
    'FontWeight','bold');
set(lh,'FontName',fontname,'FontSize',legendfontsize);
fh = gcf;
set(fh, 'Units','inches');
pos = get(fh, 'Position');        % [left bottom width height]
aspect = pos(4) / pos(3);         % height / width
new_width  = width_in;
new_height = width_in * aspect;
set(fh, 'Position', [pos(1) pos(2) new_width new_height]);


%========================================================================|
% Comparison with Trapezoidal Direct Quadrature (Structure-Preservation) |
%========================================================================|
AoI_problem.T = 8;
AoI_problem.h = 2.5e-1;
AoI_problem.beta = 6e-5;

%======================================|
% NSFD and Trapezoidal DQ computations |
%======================================|
[t, Y_test, ~] = NSFD_AoI(AoI_problem);
[t_trap,Y_trap] = Trap_AoI(AoI_problem.S0,AoI_problem.T,AoI_problem.h,...
                           AoI_problem.A,AoI_problem.phi0,AoI_problem.beta);

figure()
subplot(1,2,1)
plot(t, Y_test(1,:),'LineWidth',2,'Color',[0.96,0.47,0.16]); hold on
plot(t_trap,Y_trap(1,:),'LineWidth',2,'Color',[0.51,0.25,0.09]);
yline(0,':')
xlh = xlabel('t');
ylh = ylabel('S(t)');
lh = legend('NSFD',...
            'TrapDQ',...
            'Location','best');
axis square tight;
axh = gca;
set(axh,'FontName',fontname,'FontSize',axisfontsize);
set([xlh,ylh],'FontName',fontname,'FontSize',labelfontsize);
set(lh,'FontName',fontname,'FontSize',legendfontsize);

width_cm = 16.51;
width_in = width_cm / 2.54;
fh = gcf;
set(fh, 'Units','inches');
pos = get(fh, 'Position');        % [left bottom width height]
aspect = pos(4) / pos(3);         % height / width
new_width  = width_in;
new_height = width_in * aspect;
set(fh, 'Position', [pos(1) pos(2) new_width new_height]);

subplot(1,2,2)
plot(t, Y_test(2,:),'g','LineWidth',2); hold on
plot(t_trap,Y_trap(2,:),'Color',[0.16,0.47,0.14],'LineWidth',2);
yline(0,':')
xlh = xlabel('t');
ylh = ylabel('φ(t)');
lh = legend('NSFD',...
            'TrapDQ',...
            'Location','best');
th = sgtitle('Numerical Solutions for h = 0.25');
axis square tight
axh = gca;
set(axh,'FontName',fontname,'FontSize',axisfontsize);
set([xlh,ylh],'FontName',fontname,'FontSize',labelfontsize);
set(th,...
    'FontName',fontname,...
    'FontSize',titlefontsize,...
    'FontWeight','bold');

fh = gcf;
set(fh, 'Units','inches');
pos = get(fh, 'Position');        % [left bottom width height]
aspect = pos(4) / pos(3);         % height / width
new_width  = width_in;
new_height = width_in * aspect;
set(fh, 'Position', [pos(1) pos(2) new_width new_height]);

fprintf('\nDone.\n\n');

%=======================|
% Trapezoidal DQ method |
%=======================|
function [t,Y] = Trap_AoI(S0,T,h,A,Phi0,beta)
N = ceil(T/h);
t = linspace(0,T,N+1);
Y = zeros(2,N+1);
phi0 = Phi0(0);
sp0  = S0*phi0;
A0   = A(0);
Y(:,1) = [S0; phi0];
opts = optimoptions('fsolve',...
                    'Display','off',...
                    'FunctionTolerance',1e-12,...
                    'StepTolerance',1e-12);
for n = 1:N
    tmp1 = 0;
    tmp2 = 0;
    if n > 1
        for j = 2:n-1
            Sj   = Y(1,j);
            Phij = Y(2,j);
            prod = Sj*Phij;
            tmp1 = tmp1 - 2*prod;
            tmp2 = tmp2 + 2*A(t(n-j))*prod;
        end
    end
    tmp = [S0; Phi0(t(n))] + ...
          (beta*h/2)*( [ -sp0         + tmp1 ;
                          A(t(n))*sp0 + tmp2 ] );
    FUN = @(X) [ X(1) + (beta*h/2)*X(1)*X(2) - tmp(1);
                 X(2) - (beta*h/2)*A0*X(1)*X(2) - tmp(2) ];
    Y(:,n+1) = fsolve(FUN,Y(:,n),opts);

end
end

% ==============================================================================
% End of AoI_Final_Size.m
% ==============================================================================
