function [t, S, phi, P, Perf] = NSFD_Renewal(varargin)
%                            NSFD Method 
% Unconditionally Positive Non Standard Finite Difference (NSFD) scheme for
% the solution of a comprehensive renewal framework which includes
% different variants of the Age-of-Infection (AoI) model
% Sᵢ'(t)=-βᵢ Sᵢ(t) Vᵢ(t)
% φᵢ(t) = φᵢ₀(t) + βᵢ ∫₀ᵗ Aᵢ(t-τ)Sᵢ(τ)Vᵢ(τ) dτ
%  P(t) = P₀(t) + ∫₀ᵗ B(t-τ) Σᵣ₌₁ᴹ cᵣφᵣ(τ) dτ
% Vᵢ(t) = Σᵣ₌₁ᴹ βᵢᵣ φᵣ(t) + αᵢP(t)
%
% Version:    1.0 (part of EPITIME software) 
% Date:       April, 2026
% Authors:    M. Pezzella, G. Zanghirati
% See:        https://arxiv.org/abs/2605.00067v1
% License:    GNU GPL v.3
% Repository: https://github.com/ghitan/EPITIME
% ==============================================================================
%                         Input (=>) and Outputs (<=)
% INPUT USAGE:
% The function supports two alternative input modes:
% (1) Structured input:
%     => NSFD_Renewal(problem_struct)
%     with fields:
%       • S0, phi0, P0, A, B, beta, BETA, alpha, c, h, T, verbosity
%       • S0       : Initial vector of susceptibles
%                                               [Default: 100*[9; 9.9; 10]]
%       • verbosity: 0 to suppress warnings, 1 to enable       [Default: 1]
%
% (2) Positional input:
%     => NSFD_Renewal(S0, phi0, P0, A, B, beta, BETA, alpha, c, h, T, ...
%                                                                verbosity)
%     Missing fields are replaced by default values.
%
% OUTPUTS:
% <= t : Time vector (uniform grid with step h, except possibly for the
%        last point if T is not an integer multiple of h)
% <= S    : Susceptibles matrix (M x Nt)
% <= phi  : Infectivity matrix (M x Nt)
% <= P    : P vector (1 x Nt)
% <= Perf : (Optional) Structure with performance metrics:
%          • elapsed_time: Total execution time (seconds)
%          • steps_number: Total number of time steps performed
%          • flop_number : Estimated total number of floating-point ops.
%          • Kern_number : Number of kernel functions (A, B) evaluations
%          • memory_bytes: Estimated memory used for core variables (bytes)
% ==========================================================================
%
% For further details see:
% E. Messina, M. Pezzella, A. Vecchio.
% "Nonlocal finite difference discretization of a class of renewal equation
%  models for epidemics"
% Mathematical Biosciences and Engineering, 2023, 20(7): 11656-11675.
% DOI: 10.3934/mbe.2023518
% 
% ==========================================================================
% 
% DISCLAIMER
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% 
% ==========================================================================

%==============================|
% Input parsing and setup      |
%==============================|
Renewal_Problem = parse_input(varargin{:});

%=========================|
% Time vector setup       |
%=========================|
h = Renewal_Problem.h;  T = Renewal_Problem.T;
t = 0:h:T;
flag_small = 0;
if t(end)~=T
    % Stpsize reduction if T/h is not integer
    flag_small = 1;           h_small = T-t(end);
    t(end+1) = T;
end
Nt = length(t);
M = length(Renewal_Problem.S0);

%=========================|
% Preallocations          |
%=========================|
tic
S   = Renewal_Problem.S0 .* ones(M, Nt);
phi = zeros(M, Nt);
P   = zeros(1, Nt);
V   = zeros(M, Nt);
A_vals    = zeros(M, Nt);
phi0_vals = zeros(M, Nt);
for i = 1:M
    A_vals(i, :)    = Renewal_Problem.A(i, t);
    phi0_vals(i, :) = Renewal_Problem.phi0(i, t);
end
B_vals   = Renewal_Problem.B(t);
P0_vals  = Renewal_Problem.P0(t);
alpha    = Renewal_Problem.alpha;                BETA = Renewal_Problem.BETA;
c        = Renewal_Problem.c;                    beta = Renewal_Problem.beta;
phi(:,1) = phi0_vals(:,1);
P(1)     = Renewal_Problem.P0(t(1));
V(:,1)   = alpha .* P(1) + BETA * phi(:,1);
H        = h*beta;
allocation_time = toc;

%=========================|
% Consistency check       |
%=========================|
inputcheck_Renewal(Renewal_Problem, A_vals, B_vals, P0_vals);

%=========================|
% Time-stepping loop      |
%=========================|
tic;  % Start timing
for n = 1:Nt-1
    if  (flag_small == 1) && (n == Nt-1)
        H = Renewal_Problem.beta * h_small;
    end
    % Update of S
    S(:,n+1) = S(:,n) ./ (1 + H.*V(:,n));
    % Update of φ
    conv_phi = sum(S(:, 2:n+1) .* A_vals(:, n:-1:1) .* V(:, 1:n), 2);
    phi(:,n+1) = phi0_vals(:,n+1) + H.*conv_phi;
    % Update of P
    B_vec = B_vals(n:-1:1);
    P(n+1) = P0_vals(n+1) + h*(B_vec* (c' * phi(:, 1:n) )' );
    % Update of V
    V(:,n+1) = alpha.*P(n+1) + BETA*phi(:,n+1);
end
stepping_time = toc;
elapsed_time  = allocation_time + stepping_time;

%=======================================|
% Performance metrics (optional output) |
%=======================================|
if nargout > 4
    flop_counter = 0;
    for n = 1:Nt-1
        flop_counter = flop_counter + M*3;                    % Update of S
        flop_counter = flop_counter + M*(3*n) + 2*M;          % Update of φ
        flop_counter = flop_counter + 2*M*n + 2*n + 1;        % Update of P
        flop_counter = flop_counter + 2*M*M +2*M;             % Update of V
    end
    % Estimate memory usage (approximate core variables)
    mem_S = whos('S');
    mem_phi = whos('phi');
    mem_P = whos('P');
    mem_V = whos('V');
    mem_A_vals = whos('A_vals');
    mem_phi0_vals = whos('phi0_vals');
    mem_B_vals = whos('B_vals');
    total_mem_bytes = mem_S.bytes + mem_phi.bytes + mem_P.bytes + ...
        mem_V.bytes +  mem_A_vals.bytes + mem_phi0_vals.bytes + ...
        mem_B_vals.bytes;
    % Store performance info in output structure
    Perf = struct();
    Perf.steps_number = Nt-1;
    Perf.flop_number  = flop_counter;
    Perf.memory_bytes = total_mem_bytes;
    Perf.Kern_number  = (M+1)*Nt;
    Perf.elapsed_time = elapsed_time;
end
end
%==========================================================================
%==========================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY MODULE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% PARSE_INPUT  Parse and complete the Renewal problem input %%%%%%%%
function Renewal_Problem = parse_input(varargin)
% Handles input parsing, default assignment and verbosity setup for the
% NSFD_Renewal solver. Returns a fully initialized struct Renewal_Problem.
% Usage:
%   Renewal_Problem = parse_input(problem_struct)
%   Renewal_Problem = parse_input(S0, phi0, P0, A, B, ...)
%==============================|
% Default parameter definition |
%==============================|
defaults = struct('S0', [900; 990; 1000], 'phi0',[], 'P0', @(t) exp(-2*t), ...
    'A',@(i,t) (i*0.5) * gampdf(t,3,2), 'B',@(t) exp(-t), ...
    'beta',1e-3, 'BETA',ones(3,3) , 'alpha',[3; 1; 2],...
    'c', [0.5; 0.8; 0.7], 'h',1e-2, 'T', 1e+2, 'verbosity', 1);
N=[1000; 1050; 1010];
defaults.phi0 = @(i,t) (N(i)-defaults.S0(i))*defaults.A(i,t);
fields = fieldnames(defaults);

%======================================================|
% Parse input: either a struct or positional arguments |
%======================================================|
args = varargin;
if nargin == 1 && isstruct(args{1})
    Renewal_Problem = args{1};
    for k = 1:numel(fields)
        f = fields{k};
        if ~isfield(Renewal_Problem, f)
            Renewal_Problem.(f) = defaults.(f);
        end
    end
else
    Renewal_Problem = defaults;
    n_given = min(nargin, numel(fields));
    for k = 1:n_given
        Renewal_Problem.(fields{k}) = args{k};
    end
end

%====================================|
% Verbosity and Warnings definitions |
%====================================|
verbosity = Renewal_Problem.verbosity;
if verbosity
    if nargin == 1 && isstruct(args{1})
        for k = 1:numel(fields)
            f = fields{k};
            if ~isfield(args{1}, f)
                warning('Parameter "%s" missing: default value used.', f);
            end
        end
    else
        if nargin < numel(fields)
            missing_fields = fields(nargin+1:end);
            warning('Missing inputs. Default values will be used for: %s', ...
                strjoin(missing_fields, ', '));
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY MODULE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUTCHECK_RENEWAL Validate Renewal_Problem structure consistency %%%%
function inputcheck_Renewal(Renewal_Problem, A_vals, B_vals, P0_vals)
% It performs validation and consistency checks on scalar, vector and
% matrix parameters of the AOI problem definition
% INPUT:
%   => Renewal_Problem : struct containing problem definition
%   => A_vals          : matrix of A_i(t_n) values (M x Nt)
%   => B_vals          : vector of B(t_n) values (1 x Nt)
%   => P0_vals         : vector of P0(t_n) values (1 x Nt)
%
% CHECKS PERFORMED:
%   (1) S0    : vector with all positive entries
%   (2) beta  : positive scalar
%   (3) BETA  : matrix with non-negative entries
%   (4) alpha : vector with non-negative entries
%   (5) h, T  : positive scalars
%   (6) A_vals, B_vals, P0_vals : non-negative arrays
%
% NOTES:
%   In case of inconsistencies, execution is stopped with a descriptive
%   error message. Warnings are issued only for mild inconsistencies.

%=====================|
% Scalar fields check |
%=====================|
posScalar = @(x,fname) validateattributes(x, {'numeric'}, ...
    {'scalar','positive','finite'}, mfilename, fname);
% Positive scalar checks
posScalar(Renewal_Problem.h, 'h');
posScalar(Renewal_Problem.T, 'T');

%==========================|
% Vector and matrix checks |
%==========================|
% --- S0 ---
validateattributes(Renewal_Problem.S0, {'numeric'}, ...
    {'real','vector','nonempty','finite'}, mfilename, 'S0');
if any(Renewal_Problem.S0 <= 0)
    error('Invalid input: all components of S0 must be strictly positive.');
end
% --- beta ---
validateattributes(Renewal_Problem.beta, {'numeric'}, ...
    {'real','vector','nonempty','finite'}, mfilename, 'beta');
if any(Renewal_Problem.beta < 0)
    error('Invalid input: all components of beta must be non-negative.');
end
% --- alpha ---
validateattributes(Renewal_Problem.alpha, {'numeric'}, ...
    {'real','vector','nonempty','finite'}, mfilename, 'alpha');
if any(Renewal_Problem.alpha < 0)
    error('Invalid input: components of alpha must be non-negative.');
end
% --- BETA ---
validateattributes(Renewal_Problem.BETA, {'numeric'}, ...
    {'real','2d','nonempty','finite'}, mfilename, 'BETA');
if any(Renewal_Problem.BETA(:) < 0)
    error('Invalid input: all entries of BETA must be non-negative.');
end

%========================|
% Functions consistency  |
%========================|
validateattributes(A_vals, {'numeric'}, ...
    {'real','nonnegative','finite'}, mfilename, 'A_vals');
validateattributes(B_vals, {'numeric'}, ...
    {'real','nonnegative','finite'}, mfilename, 'B_vals');
validateattributes(P0_vals, {'numeric'}, ...
    {'real','nonnegative','finite'}, mfilename, 'P0_vals');

%===============================|
% Dimension consistency checks  |
%===============================|
M = numel(Renewal_Problem.S0);
if size(A_vals,1) ~= M
    error('Inconsistent size: A_vals must have M rows (same as length(S0)).');
end
if numel(Renewal_Problem.alpha) ~= M
    error('Inconsistent size: alpha must have length equal to M = length(S0).');
end
if size(Renewal_Problem.BETA,1) ~= M || size(Renewal_Problem.BETA,2) ~= M
    error('Inconsistent size: BETA must be an M-by-M matrix.');
end
end
