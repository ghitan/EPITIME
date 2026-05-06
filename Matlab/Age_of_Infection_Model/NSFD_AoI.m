function [t, Y, P] = NSFD_AoI(varargin)
%                           NSFD Method
% Unconditionally Positive, Monotonicity-preserving Non Standard Finite
% Difference (NSFD) scheme for the solution of the Age of Infection (AoI)
% epidemic model
% S'(t) = β S(t)φ(t)                        || S(t): number of susceptibles
% φ(t) = φ₀(t) + β ∫₀ᵗ A(t-τ) S(τ)φ(τ) dτ   || φ(t): total infectivity
% 
% Version:    1.0
% Date:       April, 2026
% Authors:    M. Pezzella, G. Zanghirati
% See:        https://arxiv.org/abs/2605.00067v1
% License:    GNU GPL v.3
% Repository: https://github.com/ghitan/EPITIME
% ==========================================================================
%                         Input (=>) and Outputs (<=)
% INPUT USAGE:
% The function supports two alternative input modes:
% 
% (1) Structured input:
%     => NSFD_AoI(AoI_problem)
%     where AoI_problem is a structure with the following optional fields:
%       • S0       : Initial number of susceptibles S(0)      [Default: 990]
%       • N        : Population Size                         [Default: 1000]
%       • T        : Final time T (integration over [0,T])     [Default: 10]
%       • h        : Time step size                           [Default: 0.1]
%       • beta     : Effective contact rate                [Default: 0.0025]
%       • A        : Infectivity kernel A(t)             [Default: exp(-3t)]
%       • phi0     : Initial infectivity φ₀(t)        [Default: (N–S0)·A(t)]
%       • verbosity: 0 to suppress warnings, 1 to enable        [Default: 1]
% 
% (2) Positional input:
%     => NSFD_AoI(S0, N, T, h, beta, A, phi0, verbosity)
%     arguments are given in order, with missing ones replaced by defaults.
%     The last parameter 'verbosity' is optional.
% 
% OUTPUTS:
% <= t : Time vector on a uniform grid with spacing h, with the last point
%        adjusted if T is not an integer multiple of h
% <= Y : Solution array. Each column is [Sₙ; φₙ] ≈ [S(tₙ); φ(tₙ)]
% <= P : (Optional) Structure with performance metrics:
%         • elapsed_time: Total execution time (seconds) excluding checks
%         • steps_number: Total number of time steps performed
%         • flops_number: Estimated number of floating-point operations
%                         not including kernel function evaluations,
%                         for which the computational cost is not assessed
%         • Kern_number : Number of kernel function evaluations
%         • memory_bytes: Estimated memory used by core variables (bytes)
% ==========================================================================
% 
% For further details see:
% E. Messina, M. Pezzella, A. Vecchio.
% "A non-standard numerical scheme for an age-of-infection epidemic model"
% Journal of Computational Dynamics, 2022, 9(2): 239–252.
% DOI: 10.3934/jcd.2021029
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

%==========================================================================
% The MATLAB implementation is organized into logical Matlab Code Blocks,
% labelled MCB:A to MCB:F, each corresponding to a specific stage of the
% solver. These labels are used consistently in the accompanying paper
% to reference individual parts of the implementation.
%
% MCB:A : Default parameters and input processing
% MCB:B : Time grid construction and memory allocation
% MCB:C : Input consistency and admissibility checks
% MCB:D : NSFD time-stepping scheme
% MCB:E : Performance metrics and diagnostics
% MCB:F : Auxiliary routine 1
% MCB:G : Auxiliary routine 2
%==========================================================================

%==========================================================================
% MCB:A | Default parameters and input processing
%--------------------------------------------------------------------------
% Definition of default values for all model parameters. These defaults
% are used whenever the corresponding input is not provided by the user,
% either in structured or positional input mode.
%==========================================================================

%==============================|
% Default parameter definition |
%==============================|
defaults = struct('S0',990, 'N',10^3, 'T',10, 'h',0.1, 'beta',2.5e-3, ...
    'A',@(t)( exp(-3*t) ), 'phi0',[], 'verbosity',1);
defaults.phi0 = @(t) ((defaults.N - defaults.S0) * defaults.A(t));
%======================================================|
% Parse input: either a struct or positional arguments |
%       Verbosity and Warnings definitions             |
%======================================================|
AoI_problem = parse_input(defaults, varargin{:});

%==========================================================================
% MCB:B | Time grid construction and memory allocation
%--------------------------------------------------------------------------
% Construction of the temporal grid, handling the case where T/h is not
% an integer. Memory is allocated once to guarantee efficiency and to
% ensure reproducibility of performance metrics.
%==========================================================================

%=====================================|
% Initialize time and solution arrays |
%=====================================|
tic;
h = AoI_problem.h;
t = 0:h:AoI_problem.T;
flag_small = 0;
if ~isempty(t)
    if t(end) ~= AoI_problem.T
        % Stpsize reduction if T/h is not integer
        flag_small = 1;
        h_small = AoI_problem.T - t(end);
        t(end+1) = AoI_problem.T;
    end
    Nt = length(t);
    H = AoI_problem.beta * h;
    P0 = AoI_problem.phi0(t);
    A_val = AoI_problem.A(t');
    % Allocate memory for solution matrix Y
    % Y(1,:) stores S(t_n), Y(2,:) stores φ(t_n)
    Y = zeros(2, Nt);
    Y(:,1) = [AoI_problem.S0; P0(1)];
else
    A_val = []; P0 = [];
end
allocation_time = toc;

%==========================================================================
% MCB:C | Input consistency and admissibility checks
%--------------------------------------------------------------------------
% Verification of positivity, finiteness and logical consistency of all
% inputs, including kernel non-negativity and admissibility of φ₀.
%==========================================================================

%=========================|
% Input Consistency Check |
%=========================|
check_input(AoI_problem, defaults, A_val, P0, AoI_problem.verbosity);

%==========================================================================
% MCB:D | NSFD time-stepping scheme
%--------------------------------------------------------------------------
% Core numerical loop implementing the NSFD discretization of the AoI
% model. The scheme guarantees unconditional positivity and monotonicity
% preservation by construction.
%==========================================================================

%=====================|
% Time-stepping loop  |
%=====================|
tic;
for n = 1:Nt-1
    if flag_small == 1 && n == Nt-1
        H = AoI_problem.beta * h_small;
    end
    % NSFD update for susceptibles
    Y(1,n+1) = Y(1,n) / (1 + H * Y(2,n));
    % Convolution term: discrete approximation of ∫₀ᵗ A(t-τ)S(τ)φ(τ) dτ
    tmp = (Y(2,1:n) .* Y(1,2:n+1)) * A_val(n+1:-1:2);
    Y(2,n+1) = P0(n+1) + H * tmp;
end
stepping_time = toc;
elapsed_time = allocation_time + stepping_time;

%==========================================================================
% MCB:E | Performance metrics and diagnostics
%--------------------------------------------------------------------------
% Optional block activated only if a third output argument is requested.
% Provides timing, memory usage and a FLOP count estimate.
%==========================================================================

%=======================================|
% Performance metrics (optional output) |
%=======================================|
if nargout > 2
    % Approximate flops count:
    % 1 multiplication for H
    % Each iteration does:
    %  1 sum, 1 multiplication, 1 division for S update
    %  n multiplications + n multiplications and n-1 additions for tmp
    %  2 additions and 1 multiplication for φ update
    % Total ≈ 6+3n-1 flops per step
    flop_counter = 1 + sum( 5 +3 * (1:Nt-1) ) + 2*flag_small;
    % Estimate memory usage of core variables
    mem_Y = whos('Y');
    mem_t = whos('t');
    mem_A_val = whos('A_val');
    mem_P0 = whos('P0');
    total_memory_bytes = mem_Y.bytes + mem_t.bytes + mem_A_val.bytes ...
        + mem_P0.bytes;
    % Store performance info in output structure
    P = struct();
    P.elapsed_time = elapsed_time;
    P.steps_number = Nt - 1;
    P.flop_number = flop_counter;
    P.Kern_number = Nt;
    P.memory_bytes = total_memory_bytes;
end
end
%==========================================================================
%==========================================================================

%==========================================================================
% MCB:F and MCB:G | Auxiliary routines
%--------------------------------------------------------------------------
% The following local functions implement modular tasks that are logically
% separated from the main solver, namely input parsing and validation.
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY MODULE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AoI_problem = parse_input(defaults, varargin)
%==========================================================================
% MCB:F | Parse and complete the AoI problem input
%--------------------------------------------------------------------------
% Handles input parsing, default assignment and verbosity setup for the 
% NSFD_AoI solver. Returns a fully initialized struct AoI_problem.
% Usage:
%   AoI_problem = parse_input(defaults, varargin{:})
% Handles both:
%   1. Structured input: NSFD_AoI(AoI_problem)
%   2. Positional input: NSFD_AoI(S0, N, T, h, beta, A, phi0, verbosity)
%==========================================================================
fields = fieldnames(defaults);
args = varargin;
%==========================|
% Structured input         |
%==========================|
if nargin == 2 && isstruct(args{1})   % Structured input
    AoI_problem = args{1};
    missing_fields = {};
    % Fill in missing fields with default values
    for k = 1:numel(fields)
        f = fields{k};
        if ~isfield(AoI_problem, f)
            AoI_problem.(f) = defaults.(f);
            missing_fields{end+1} = f; 
        end
    end
    % Emit a single warning if needed
    verbosity = AoI_problem.verbosity;
    if verbosity && ~isempty(missing_fields)
        warning('Missing inputs. Default values will be used for: %s.', ...
                strjoin(missing_fields, ', '));
    end
    %==========================|
    % Positional input         |
    %==========================|
else
    AoI_problem = defaults;
    n_given = min(length(args), numel(fields));
    for k = 1:n_given
        if ~isempty(args{k})
            AoI_problem.(fields{k}) = args{k};
        end
    end
    % Generate default phi0 if empty
    if isempty(AoI_problem.phi0)
        AoI_problem.phi0 = @(t) (AoI_problem.N - AoI_problem.S0) ...
            * AoI_problem.A(t);
    end
    % Warnings for missing positional arguments
    verbosity = AoI_problem.verbosity;
     if verbosity
        missing_fields = {};
        for k = 1:numel(fields)
            if isequal(AoI_problem.(fields{k}), defaults.(fields{k})) ...
               && (k > length(args) || isempty(args{k}))
                missing_fields{end+1} = fields{k}; 
            end
        end
        if ~isempty(missing_fields)
            warning('Missing inputs. Default values will be used for: %s.',...
                    strjoin(missing_fields, ', '));
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY MODULE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_input(AoI_problem, defaults, A_val, P0, verbosity)
%==========================================================================
% MCB:G | Validate AoI_problem structure fields and consistency
%--------------------------------------------------------------------------
% It performs validation and consistency checks on scalar, vector and
% matrix parameters of the AoI problem definition
% INPUT:
%   => AoI_problem : struct containing problem definition
%   => defaults    : struct with default values for each field
%   => A_val       : vector of A(t_n) values (for phi0 check)
%   => P0          : vector phi0(t_n) initial values
%   => verbosity   : 0 to suppress warnings, 1 to enable
% CHECKS PERFORMED:
%   (1) Positivity and finiteness of scalar fields
%   (2) Non-negativity of A(t) and phi0(t)
%   (3) Logical relation: N > S0
%   (4) Elementwise consistency: phi0(t_n) ≤ (N - S0)*A(t_n)
%
% NOTES:
%   Errors stop execution. Warnings are issued only for mild inconsistencies
%   (phi0 exceeding its theoretical upper bound).
%==========================================================================

%=====================|
% Scalar fields check |
%=====================|
% Anonymous validator for positive scalars
posScalar = @(x,fname) validateattributes(x, {'numeric'}, ...
    {'scalar','positive','finite'}, ...
    mfilename, fname);
scalar_fields = {'S0', 'N', 'T', 'h', 'beta'};
% Check that numerical inputs (if user-defined) are positive and consistent
for k = 1:numel(scalar_fields)
    fname = scalar_fields{k};
    value = AoI_problem.(fname);
    user_defined = isfield(AoI_problem, fname) && ...
        ~isequal(value, defaults.(fname));
    if user_defined
        posScalar(value, fname);
    end
end
%========================================================================|
% Functions checks (Non-negativity, logical relations and φ₀ upper bound |
%========================================================================|
validateattributes(A_val, {'numeric'}, {'real','vector','nonnegative'}, ...
    mfilename, 'A_val');
validateattributes(P0, {'numeric'}, {'real','vector','nonnegative'}, ...
    mfilename, 'phi0');
% Ensure N > S0
if AoI_problem.N <= AoI_problem.S0
    error('Invalid values: N must be greater than S0.');
end
if isfield(AoI_problem, 'phi0') && ...
        ~isequal(AoI_problem.phi0, defaults.phi0)
    validateattributes(P0, {'numeric'}, {'vector','real','nonnegative'}, ...
        mfilename, 'phi0');
    upper_bound = (AoI_problem.N - AoI_problem.S0) * A_val(:);
    % Check elementwise consistency φ₀(tₙ) ≤ (N-S0)A(tₙ)
    bound_sensitivity = 2e-14;
    exceeding = P0(:) > upper_bound + bound_sensitivity;
    if any(exceeding)
        if verbosity
            warning(['Inconsistent phi0: some values exceed (N - S0)*A(t).' ...
                ' Check your definition of phi0.']);
        end
    end
end
end
