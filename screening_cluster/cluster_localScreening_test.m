taskID = 2;
ID = ceil(taskID / 7);
nP = taskID - (ID - 1) * 7;
foldername = '\screen_local2\';
filename = [pwd foldername 'P' num2str(nP) '_' num2str(ID, '%05d') '.mat'];
njobs = 2000;

if ~exist(filename, 'file')

% reseed the random number generator
rng('shuffle'); rngState = rng;
seed = rngState.Seed + uint32(feature('getpid')); rng(seed);

pyramid_initRatios = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; 0.9 0.1 0; 0.99 0.01 0; 0.1 0.9 0];
Output_Biomass = zeros(3 * length(pyramid_initRatios), 1);
Output_Sizes   = zeros(3 * length(pyramid_initRatios), 1);

perturbmag = 0.2; % magnitude of perturbation to parameters
load('selectedParameters_local1.mat','ParameterSelected')
Parameters_origin = ParameterSelected(:, nP);
Parameters = Parameters_origin .* (1 + perturbmag *(2 * rand(size(ParameterSelected,1), 1) - 1)); % perturb each parameter p within the [1-perturbmag, 1+perturbmag]*p range
Parameters = roundn(Parameters, -1);
if ID == 1; Parameters = Parameters_origin; end

% Parameters
L      = 90;    % domain size
totalt = 16;    % total time
dt     = 0.02;  % time step
nx     = 1001; ny = nx; % number of nodes

bN = Parameters(1);         % nutrient consumption rate
DN = Parameters(2);          % nutrient diffusivity
N0 = Parameters(3);         % initial nutrient conc.

aCs_act = [1, Parameters(7), 1] * Parameters(4);  % cell growth rate of each species
gs  = [1, 1, Parameters(8)] * Parameters(5);        % swimming motility
hs_act  = [1, 0, 0.9] * Parameters(6);  % swarming motility coefficients
mu = Parameters(9);
N_upper = Parameters(10) * N0; % upper bound of nutrient for swarming
N_lower = N_upper - Parameters(11) * N0; % lower bound of nutrient for swarming

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BranchingColonyMultispecies_core
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save(filename, 'Parameters','Output_Biomass','Output_Sizes')

end

% Collect data
if mod(taskID, 1) == 0
ParametersCell = cell(size(ParameterSelected, 2), 1);
OutputCell_Biomass = cell(size(ParameterSelected, 2), 1);
OutputCell_Sizes = cell(size(ParameterSelected, 2), 1);
ncompleted = 0;
for nP = 1 : size(ParameterSelected, 2)    
    ParametersMat = nan(length(Parameters), njobs);
    OutputMat_Biomass = nan(3 * length(pyramid_initRatios), njobs);
    OutputMat_Sizes = nan(3 * length(pyramid_initRatios), njobs);
    for id = 1 : njobs
        filename = [pwd foldername 'P' num2str(nP) '_' num2str(id, '%05d') '.mat'];
        if exist(filename, 'file')
            load(filename)  
            ParametersMat(:, id) = Parameters;
            OutputMat_Biomass(:, id) = Output_Biomass;
            OutputMat_Sizes(:, id) = Output_Sizes;
        end
    end
    ParametersCell{nP} = ParametersMat;
    OutputCell_Biomass{nP} = OutputMat_Biomass;
    OutputCell_Sizes{nP} = OutputMat_Sizes;
    ncompleted = ncompleted + sum(~isnan(ParametersMat(1,:)));
end
save('output_local_2.mat','ncompleted','ParametersCell','OutputCell_Biomass','OutputCell_Sizes')
end