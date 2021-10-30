taskID = str2double(getenv('SLURM_ARRAY_TASK_ID'));
foldername = '/screen_coarse/';
filename = [pwd foldername num2str(taskID, '%05d') '.mat'];
njobs = 20000;

if ~exist(filename, 'file')

% reseed the random number generator
rng('shuffle'); rngState = rng;
seed = rngState.Seed + uint32(feature('getpid')); rng(seed);

pyramid_initRatios = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; 0.9 0.1 0; 0.99 0.01 0; 0.1 0.9 0];
Output_Biomass = zeros(3 * length(pyramid_initRatios), 1);
Output_Sizes   = zeros(3 * length(pyramid_initRatios), 1);

% Parameters
L      = 90;    % domain size
totalt = 16;    % total time
dt     = 0.02;  % time step
nx     = 1001; ny = nx; % number of nodes

bN   = rand * (50 - 10) + 10;  % nutrient consumption rate
DN   = rand * (25 - 5) + 5;  % nutrient diffusivity
N0   = rand * (50 - 5) + 5;  % initial nutrient conc.
aCs0 = rand * (1.5 - 0.5) + 0.5;
gs0  = rand * (5 - 0.5) + 0.5;
hs0  = rand * (10 - 1) + 1;
aC2  = rand * (1.5 - 1.05) + 1.05;
gs3  = rand * (3 - 1.5) + 1.5;
mu   = rand * (1 - 0.2) + 0.2;
wint = rand * (1 - 0.4) + 0.4;
winh = rand * (0.4 - 0.1) + 0.1;

aCs_act = [1, aC2, 1] * aCs0;  % cell growth rate of each species
gs  = [1, 1, gs3] * gs0;        % swimming motility
hs_act  = [1, 0, 0.9] * hs0;  % swarming motility coefficients
N_upper = wint * N0; % upper bound of nutrient for swarming
N_lower = N_upper - winh * N0; % lower bound of nutrient for swarming

Parameters = [bN, DN, N0, aCs0, gs0, hs0, aC2, gs3, mu, wint, winh];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BranchingColonyMultispecies_core
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save(filename, 'Parameters','Output_Biomass','Output_Sizes')

end

% Collect data
if mod(taskID, 100) == 0  
    % collect data
    ParametersMat = nan(length(Parameters), njobs);
    OutputMat_Biomass = nan(3 * length(pyramid_initRatios), njobs);
    OutputMat_Sizes = nan(3 * length(pyramid_initRatios), njobs);
    for id = 1 : njobs
        filename = [pwd foldername num2str(id, '%05d') '.mat'];
        if exist(filename, 'file')
            load(filename)  
            ParametersMat(:, id) = Parameters;
            OutputMat_Biomass(:, id) = Output_Biomass;
            OutputMat_Sizes(:, id) = Output_Sizes;
        end
    end
    ncompleted = sum(~isnan(ParametersMat(1,:)));
    save('output.mat','ncompleted','ParametersMat','OutputMat_Biomass','OutputMat_Sizes')
end