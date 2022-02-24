speciesName = {'WT','Cheater','Hyperswarmer'}; % name of each species
% other vectors will follow the same order

%% Tunable parameters

% bN = 25;         % nutrient consumption rate
% DN = 7;          % nutrient diffusivity
% N0 = 20;         % initial nutrient conc.
% 
% aCs_act = [1, 1.1, 1] * 1.5/2;      % cell growth rate of each species
% gs  = [1, 1, 2] * 1;            % swimming motility
% hs_act  = [1, 0, 0.9] * 8;        % swarming motility coefficients
% 
% N_upper = 15; % upper bound of nutrient for swarming
% N_lower = 6; % lower bound of nutrient for swarming

%% Imported parameters
if length(Parameters) < 12; Parameters(12) = 1; end
if length(Parameters) < 13; Parameters(13) = 0.9; end
bN = Parameters(1);         % nutrient consumption rate
DN = Parameters(2);          % nutrient diffusivity
N0 = Parameters(3);         % initial nutrient conc.

aCs_act = [1, Parameters(7), 1] * Parameters(4);  % cell growth rate of each species
gs  = [1, 1, Parameters(8)] * Parameters(5);        % swimming motility
hs_act  = [1, 0, 0.9] * Parameters(6);  % swarming motility coefficients

N_upper = Parameters(10) * N0; % upper bound of nutrient for swarming
N_lower = N_upper - Parameters(11) * N0; % lower bound of nutrient for swarming

%% Fixed parameters
L      = 90;    % domain size
totalt = 16;    % total time
dt     = 0.02;  % time step
nx     = 1001; ny = nx; % number of nodes

r0  = 5;   % initial radius
C0  = 8;   % initial cell density

% branch density & width of single-species colonies
Densities = [0.14, 0.14, 0.25];
Widths    = [5, 5, 20];

noiseamp = 0;        % noise amplitude of branch direction
dt_updatebranch = 10 * dt;  % time step for updating branch locations
