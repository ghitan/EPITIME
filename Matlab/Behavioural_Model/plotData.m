% Plot of Test1_JMB for the paper experiments (part of EPITIME software)
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

%============================================|
% Plot of the numerical solution  for h=1e-1 |
%============================================|
function plotData(t,F,Feq)

% -----------------------------
% Impostazioni tipografiche
% -----------------------------
axisfontsize   = 11;
labelfontsize  = 11;
titlefontsize  = 13;
legendfontsize = labelfontsize;

%fontname = 'Palatino Linotype';   % Windows
 fontname = 'Palatino';          % macOS
% fontname = 'URW Palladio L';     % Linux

% -----------------------------
% Figura
% -----------------------------
figure
%plot(t, F./Feq, 'b', 'LineWidth', 2); hold on
plot(t, F./Feq, 'Color', [0 0.4 0], 'LineWidth', 2); hold on

% Linea orizzontale y = 1 (azzurra tratteggiata)
%yline(1, '--', 'Color', [0 0.6 0.8], 'LineWidth', 1.5);
yline(1, '--', 'Color', [0.4 0.8 0.4], 'LineWidth', 1.5);


xlh = xlabel('t');
ylh = ylabel('$$F(t)\,/\,F_e$$','Interpreter','latex');
th  = title('NSFD Numerical Solution');

lh = legend({'$$F/F_e$$', '$$F/F_e=1$$'}, ...
            'Interpreter','latex', ...
            'Location','northeast');

axis tight;

axh = gca;

% -----------------------------
% Stile assi
% -----------------------------
set(axh, ...
    'FontName', fontname, ...
    'FontSize', axisfontsize);

set([xlh, ylh], ...
    'FontName', fontname, ...
    'FontSize', labelfontsize);

set(th, ...
    'FontName', fontname, ...
    'FontSize', titlefontsize);

set(lh, ...
    'FontName', fontname, ...
    'FontSize', legendfontsize);

% -----------------------------
% Dimensione fisica (full width)
% -----------------------------
width_cm = 16.51;
width_in = width_cm / 2.54;

fh = gcf;
set(fh, 'Units','inches');

pos = get(fh, 'Position');        % [left bottom width height]
aspect = pos(4) / pos(3);         % height / width

new_width  = width_in;
new_height = width_in * aspect;

set(fh, 'Position', [pos(1) pos(2) new_width new_height]);

% -----------------------------
% Export vettoriale
% -----------------------------
fname = 'Test_JMB_beh_OSCI_trap.pdf';
exportgraphics(fh, fname, ...
    'BackgroundColor','none', ...
    'Resolution',1200, ...
    'ContentType','vector');
end
