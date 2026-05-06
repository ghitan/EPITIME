function [t, Y, P] = NSFD_behavioural(varargin)
%                                   NSFD_behavioural
% ____________________________________________________________________________________________
%   Weighted Non Standard Finite Difference (NSFD) scheme for the solution of the behavioural
%   integro-differential epidemic model with demographic turnover and information index
%
% Version:    1.0 (part of EPITIME software)
% Date:       April, 2026
% Authors:    C. Panico, G. Zanghirati
% See:        https://arxiv.org/abs/2605.00067v1
% License:    GNU GPL v.3
% Repository: https://github.com/ghitan/EPITIME
% ********************************************************************************************
%   SYNOPSIS:
%   [t, Y, P] = NSFD_behavioural(varargin)
%
% ********************************************************************************************
%   INPUT:
%   The function supports two alternative input modes:
%
%   1) The function has a structured input BEH_prob with the following optional fields:
%       
%       T             : Final time of the integration                    [Default: 1000]
%       h             : Time step used for discretization                [Default: 1]
%       mu            : Natural death rate                               [Default: 1/(75*365)]
%       N             : Total population size                            [Default: 5e7]
%       S0            : Initial susceptible population                   [Default: 0.2*6e7]
%       A             : Infectivity function                             [Default: bb* t* exp(-(mu +ni).*t)]
%       betaM         : Inhibition function                              [Default: 1/((1 / N) + alpha* x)]
%       g             : Message function                                 [Default: x)]
%       K             : Memory kernel                                    [Default: a * exp(-a * x)]
%       Verbosity     : 0 to suppress warnings, 1 to enable              [Default: 1]
%      
%       EXAMPLE: NSFD_behavioural(BEH_prob)
%     
%     2) The function has positional inputs:
%      
%       EXAMPLE: NSFD_behavioural(T, h, N, mu, S0, A, K, betaM, g, Verbosity)
%
%  OUTPUT:
%   t                   : vector of time 
%   s =Y(1,:)           : numerical approximation of the solution s(t)
%   f =Y(2,:)           : numerical approximation of the solution f(t) 
%   P                   : (optional) Struct containing performance metrics and resource usage:
%                          • elapsed_time   - total computation time [seconds]
%                          • steps_number   - number of time steps or iterations (Nt)
%                          • flops_number   - estimated number of floating-point operations
%                          • memory_bytes   - total memory used by key variables [bytes]
%
% _________________________________________________________________________________________
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
% _________________________________________________________________________________________
%
%  For further details we refer to the following papers
%  [1]  B. Buonomo, E. Messina, C. Panico, A. Vecchio.
%        A stable numerical method for integral epidemic models with behavioral changes in 
%        contact patterns. 
%        Electronic Transactions on Numerical Analysis, 2024, 61, 137–156.
%        DOI : 10.1553/etna_vol61s137
%
%  [2]  B. Buonomo, E. Messina, C. Panico, A. Vecchio.
%        An integral renewal equation approach to behavioural epidemic models with 
%        information index
%        Journal of Mathematical Biology, 2025, 90(8).
%        DOI : 10.1007/s00285-024-02172-y
% _________________________________________________________________________________________
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
% _________________________________________________________________________________________

% _________________________________________________________________________________________
% MCB:A | Default parameters and input processing
%
% ----- Default parameters -----
defaults = struct( ...
    'T',         1000, ...
    'h',         1, ...
    'N',         5e7, ...
    'mu',        1/(75*365), ...
    'S0',        0.2*5e7, ...
    'A',         [], ...
    'K',         [], ...
    'betaM',     [], ...
    'g',         [], ...
    'verbosity', 1 ...
    );
% ----- Input parsing new -----
BEH_prob = parse_input(defaults, varargin{:});
% _________________________________________________________________________________________

% _________________________________________________________________________________________
% MCB:B | Time grid construction and memory allocation
%
% ---------- Model Setup Starts Here ------------------
tic
T   = BEH_prob.T;
h   = BEH_prob.h;
N   = BEH_prob.N;
mu  = BEH_prob.mu;
S0  = BEH_prob.S0 / N;


ni_t    = floor(T / h);
D       = h;
D_tilde = 0;
t       = 0:h:T;

if ( t(end) ~= T )
   D_tilde   = T-ni_t*h;
   t(ni_t+1) = T;
end
Nt = length(t);

% ------------ Memory allocation -------------------------
Y = zeros(2, Nt);
M = zeros(1, Nt);
betaM_vec = zeros(1, Nt);
g_vec     = zeros(1, Nt);

%-------------- Initial conditions -----------------------
Y(1,1)   = S0;
A_vec    = BEH_prob.A(t);
K_vec    = BEH_prob.K(t);
f_0_cost = (1 - S0/N) * mu; 
F0_vec   = f_0_cost .*A_vec;
Y(2,1)   = F0_vec(1);
M0_vec   = (1/N)*BEH_prob.g(N*Y(2,1)) .*K_vec;
M(1)     = M0_vec(1);

allocation_time = toc;
% _________________________________________________________________________________________

% _________________________________________________________________________________________
% MCB:C | Input consistency and admissibility checks
%-------Input Consistency Check---
check_input(BEH_prob, defaults, A_vec);
% _________________________________________________________________________________________

% _________________________________________________________________________________________
% MCB:D | NSFD weight definition and time-stepping scheme
%------ NSFD weight------
tic
Sum1   = sum(K_vec);
T_incr = T+h;
incr   = h*BEH_prob.K(T_incr);
it_max = 0;
Nt_2   = Nt^2;
K_eval = Nt;
while ( (incr > 1e-11) && (it_max < Nt_2) )
    Sum1   = Sum1+incr;
    T_incr = T_incr+h;
    incr   = h*BEH_prob.K(T_incr);
    it_max = it_max+1;
    K_eval = K_eval+1;
end
gammah = 1 / (h * Sum1);
weight_definition_time=toc;

% ----- Time iteration -----
tic
for n = 1:Nt-1
    if (D_tilde ~= 0 && n == Nt-1)
        D = D_tilde;
    end
    betaM_val    = BEH_prob.betaM(N*M(n));
    g_val        = BEH_prob.g(N*Y(2,n)); 
    Y(1,n+1)     = (Y(1,n) + mu*D) / (1 + D*(mu + N*betaM_val * Y(2,n)));
    betaM_vec(n) = N*betaM_val;
    g_vec(n)     = g_val;
    
    var1 = Y(1,2:n+1) .* betaM_vec(1:n) .* Y(2,1:n);
    int1 = A_vec(n+1:-1:2) * var1';
    int2 = K_vec(n:-1:1) * g_vec(1:n)';

    Y(2,n+1) = F0_vec(n+1) + D * int1;
    M(n+1)   = M0_vec(n+1) + gammah * (1/N) * D * int2;
end

stepping_time = toc;
elapsed_time  = allocation_time + weight_definition_time + stepping_time;
% _________________________________________________________________________________________

% _________________________________________________________________________________________
% MCB:E | Performance metrics and diagnostics
% ----- Optional output -----
if ( nargout > 2 )
    flops_counter = sum( 9 + 6 * (1:Nt-D) ) + 18 * D;
    mem_Y = whos('Y'); mem_t = whos('t');
    mem_A = whos('A_vec'); mem_F0 = whos('F0_vec');
    memory_bytes  = mem_Y.bytes + mem_t.bytes + mem_A.bytes + mem_F0.bytes;
    % Store performance info in output structure
    P = struct();
    P.elapsed_time = elapsed_time;
    P.steps_number = Nt;
    P.flops_number = flops_counter;
    P.Kern_number  = Nt+K_eval;
    P.memory_bytes = memory_bytes;
end
end
% _________________________________________________________________________________________


% _________________________________________________________________________________________
% MCB:F and MCB:G | Auxiliary routines
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY MODULE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PARSE_INPUT  Parse and complete the BEH problem input %%%%%%%%
function BEH_prob = parse_input(defaults, varargin)

param_names = fieldnames(defaults);
% ----- Input parsing -----
if nargin == 0
    BEH_prob = defaults;
elseif nargin == 1 && isstruct(varargin{1})
    BEH_prob = varargin{1};
else
    BEH_prob = defaults;
    n_args = min(length(varargin), length(param_names));
    for k = 1:n_args
        if ~isempty(varargin{k})
            BEH_prob.(param_names{k}) = varargin{k};
        end
    end
end

% ----- Extract verbosity -----
verbosity = BEH_prob.verbosity;

% ----- Fill in missing fields with defaults -----
for k = 1:numel(param_names)
    f = param_names{k};
    if ~isfield(BEH_prob, f) || isempty(BEH_prob.(f))
        BEH_prob.(f) = defaults.(f);
        if verbosity == 1
            warning('Parameter "%s" missing or empty: default used.', f);
        end
    end
end

% ------ Function definitions ------
ni=1/7;
R0=20;
bb = R0 * (BEH_prob.mu + ni)^2; 

if isempty(BEH_prob.A)
    BEH_prob.A = @(x) bb .* x .* exp(-(BEH_prob.mu + ni).*x);
end
if isempty(BEH_prob.K)
    a=1/30;
    BEH_prob.K = @(x) a .* exp(-a .* x);
end
if isempty(BEH_prob.betaM)
    alpha=8e3;
    BEH_prob.betaM = @(x) 1 ./ ((1 / BEH_prob.N) + alpha.* x);
end
if isempty(BEH_prob.g)
    BEH_prob.g = @(x) x;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY MODULE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CHECK_INPUT Validate BEH_prob structure fields and consistency %%%%%
function check_input(BEH_prob, defaults, A_vec)
% It performs validation and consistency checks on scalar, vector and
% matrix parameters of the BEH problem definition
% INPUT:
%   => BEH_prob          : struct containing problem definition
%   => defaults          : struct with default values for each field
%   => A_vec             : vector of A(t_n) values 
%
% CHECKS PERFORMED:
%   (1) Positivity and finiteness of scalar fields
%   (2) Non-negativity of A(t)
%   (3) Logical relation: N > S0
%
% NOTES:
%   Errors stop execution.
%------Scalar fields check 
% Anonymous validator for positive scalars
posScalar = @(x,fname) validateattributes(x, {'numeric'}, ...
    {'scalar','positive','finite'}, ...
    mfilename, fname);
scalar_fields = {'S0', 'N', 'T', 'h', 'mu'};
% Check that numerical inputs (if user-defined) are positive and consistent
for k = 1:numel(scalar_fields)
    fname = scalar_fields{k};
    value = BEH_prob.(fname);
    user_defined = isfield(BEH_prob, fname) && ...
        ~isequal(value, defaults.(fname));
    if user_defined
        posScalar(value, fname);
    end
end

% Function checks (Non-negativity)                                       

validateattributes(A_vec, {'numeric'}, {'real','vector','nonnegative'}, ...
    mfilename, 'A_vec');
% Ensure N > S0
if BEH_prob.N <= BEH_prob.S0
    error('Invalid values: N must be greater than S0.');
end
end
% _________________________________________________________________________________________