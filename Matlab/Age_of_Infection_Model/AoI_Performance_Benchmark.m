% =========================================================================
%  COMPUTATIONAL COMPLEXITY ANALYSIS OF NSFD_AoI (part of EPITIME software)
%
%  This script evaluates the computational cost of the NSFD_AoI solver by
%  measuring:
%
%   1) CPU execution time using MATLAB's TIMEIT routine.
%   2) Memory consumption of the main variables returned by the solver.
%
%  The analysis is repeated for increasing numbers of time steps Nt.
%  The resulting measurements are used to assess the empirical time and
%  memory complexity of the implementation.
%
%  Output:
%   - Semilog plots of memory and CPU time versus Nt.
%   - A summary table printed in the command window for Nt multiples of 100.
%
% ========================================================================

close all; clear; clc;
fprintf('\nEPITIME AoI: computational complexity analysis\n');

% =========================================================================
%  TEST CONFIGURATION
%
%  N_tried contains all values of Nt used in the benchmark.
%  For each Nt, the solver is executed over the fixed time interval [0,T].
% ========================================================================

N_tried = 1000:5:10000;

AoI_Problem.T = 100;
AoI_Problem.verbosity = 0;

% =========================================================================
%  PREALLOCATION
%
%  Preallocate vectors used to store the measured CPU time and memory
%  consumption. Preallocation avoids repeated memory reallocations inside
%  the benchmark loop.
% ========================================================================

executionTime   = zeros(size(N_tried));
memoryAllocated = zeros(size(N_tried));

% =========================================================================
%  WARM-UP EXECUTION
%
%  The first execution of a MATLAB function may include one-time overhead
%  due to JIT compilation, cache initialization, memory allocation, etc.
%  A preliminary call is therefore performed and discarded before starting
%  the benchmark.
% ========================================================================

AoI_Problem.h = AoI_Problem.T/5;
[~,~,~] = NSFD_AoI(AoI_Problem);

% =========================================================================
%  MAIN BENCHMARK LOOP
%
%  For each value of Nt:
%
%   1) The time step size is updated.
%   2) The execution time is measured using TIMEIT.
%   3) The solver is executed once more to retrieve memory statistics.
%
% ========================================================================

Nsize = numel(N_tried);
% Prepare the running echo to console
Ndisp = max(1., Nsize / 10.) * ( 1 : min(Nsize, 10) )';
k = 1;
fprintf('\nPerformance benchmark started...\n0%%..');

for n = 1:Nsize 
    % Current number of time steps
    nt = N_tried(n);
    % Corresponding time step size
    AoI_Problem.h = AoI_Problem.T/nt;
    % CPU TIME MEASUREMENT
    func_handle = @() NSFD_AoI(AoI_Problem);
    % Suppress possible TIMEIT warnings for extremely fast executions
    ws = warning('off','all');
    executionTime(n) = timeit(func_handle);
    warning(ws);
    % MEMORY MEASUREMENT
    % Execute the solver once and extract memory information
    [~,~,P] = NSFD_AoI(AoI_Problem);
    memoryAllocated(n) = P.memory_bytes;
    % Display running iterations percentage
    if ( n >= Ndisp(k) ), fprintf('%d%%..', fix(100. * Ndisp(k)/ Nsize)); k = k + 1; end
end
fprintf('%c%c%\nBenchmark finished.\n', char(8), char(8));

% =========================================================================
%  SAVE RESULTS
% =========================================================================
save('AOI_Benchmark_Results.mat')
%%
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

% =========================================================================
%  COMPLEXITY PLOTS
%
%  Left panel:
%     Memory consumption versus Nt.
%
%  Right panel:
%     CPU execution time versus Nt.
%
%  Reference curves are added to visually compare the measured growth
%  against theoretical complexity laws.
% ========================================================================
figure
% MEMORY COMPLEXITY
subplot(1,2,1)
semilogy(N_tried, memoryAllocated,'LineWidth',2,'Color','b');
hold on
% Reference O(Nt) curve
semilogy(N_tried,...
       0.8*N_tried.*(memoryAllocated(end)/N_tried(end)),...
       '--k');
xlh = xlabel('N_t');
ylh = ylabel('Memory (bytes)');
lh = legend('Variables Memory',...
            'O(N_t)',...
            'Location','southeast');
th = title('Memory complexity analysis');
axis square tight
axh = gca;
set(axh,'FontName',fontname,'FontSize',axisfontsize);
set([xlh,ylh],'FontName',fontname,'FontSize',labelfontsize);
set(th,'FontName',fontname,'FontSize',titlefontsize);
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
yticks(10.^(-16:1:6))
% TIME COMPLEXITY
subplot(1,2,2)
semilogy(N_tried, executionTime,'LineWidth',2,'Color','c');
hold on
% Reference O(Nt^2) curve
semilogy(N_tried,...
    0.75*(N_tried.^2).*(executionTime(end)/(N_tried(end)^2)),...
       '-.k');
xlh = xlabel('N_t');
ylh = ylabel('CPU Time (s)');
lh = legend('Measured CPU time',...
            'O(N_t^2)',...
            'Location','southeast');
th = title('Time complexity analysis');
axis square tight
axh = gca;
set(axh,'FontName',fontname,'FontSize',axisfontsize);
set([xlh,ylh],'FontName',fontname,'FontSize',labelfontsize);
set(th,'FontName',fontname,'FontSize',titlefontsize);
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
yticks(10.^(-16:1:6))

% =========================================================================
%  SUMMARY TABLE
%
%  Print benchmark results only for Nt values that are multiples of 100.
%  This provides a compact overview of the measured computational cost.
% ========================================================================
idx = mod(N_tried,100) == 0;
fprintf('\n');
fprintf('-------------------------------------------------\n');
fprintf('   Nt      CPU Time (s)    Memory (MB)\n');
fprintf('-------------------------------------------------\n');
for k = find(idx)
    fprintf('%6d    %12.4e    %12.4f\n', ...
        N_tried(k), ...
        executionTime(k), ...
        memoryAllocated(k)/(1024^2));
end
fprintf('-------------------------------------------------\n');

fprintf('\nDone.\n\n');

% ==============================================================================
% End of AoI_Performance_Benchmark.m
% ==============================================================================
