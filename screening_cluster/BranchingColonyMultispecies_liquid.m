% ID = 9351; 
% ID = 2;
Parameters = ParametersMat(:, ID);
pyramid_initRatios = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; 0.9 0.1 0; 0.99 0.01 0; 0.1 0.9 0];
WC_ratios = [0.001,0.01,0.1,10]';
WH_ratios = [0.001,0.01,0.1,10]';
pyramid_initRatios = [pyramid_initRatios; 1./(1+WC_ratios),WC_ratios./(1+WC_ratios),zeros(length(WC_ratios),1)];
pyramid_initRatios = [pyramid_initRatios; zeros(length(WH_ratios),1), WH_ratios./(1+WH_ratios),1./(1+WH_ratios)];

Output_Biomass_re = zeros(3 * length(pyramid_initRatios), 1);
Output_Sizes_re   = zeros(3 * length(pyramid_initRatios), 1);

% Parameters

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

for iter = 4 : length(pyramid_initRatios)

speciesName = {'WT','Cheater','Hyperswarmer'}; % name of each species
% other vectors will follow the same order

initialRatio = pyramid_initRatios(iter, :);   % initial ratio of all species
initialFract = initialRatio / sum(initialRatio); % initial fraction of each species

C0  = 0.001;   % initial cell density

%% Initialization

nt = totalt / dt;  % number of time steps

C = cell(3, 1); C(:) = {zeros(nt + 1, 1)};   % Cell density
N = zeros(nt + 1, 1); N(1) = N0;                   % Nutrient

for j = 1 : 3
    C{j}(1) = C0 * initialFract(j);
end

aCs = cell(3, 1);


for i = 1 : nt
    
    fN = N(i) ./ (N(i) + 1) .* (1 - (C{1}(i) + C{2}(i) + C{3}(i)));
    
    for j = 1 : 3  % j: index of species
        
        % -------------------------------------
        % Cell growth
        aCs{j} = aCs_act(j); % growth rate as a variable of time and location
        if N(i) > N_upper || N(i) < N_lower; aCs{j} = max(aCs_act); end% if nutrient is not within the range grow at full speed
        C{j}(i+1) = C{j}(i) + aCs{j} .* fN .* C{j}(i) * dt;
    
    end
    
    % -------------------------------------
    % Nutrient distribution
    dN = - bN * fN .* (aCs{1} .* C{1}(i) + aCs{2} .* C{2}(i) + aCs{3} .* C{3}(i)); % Nutrient consumption
    N(i+1)  = N(i) + dN * dt;

    % -------------------------------------
    % Cell allocation and branch extension
    

end

Output_Biomass_re((iter - 1) * 3 + 1 : iter * 3) = [sum(C{1},'all'); sum(C{2},'all'); sum(C{3},'all')];
fprintf('iter = %d\n', iter)

save([picname '_' num2str(iter,'%02d') '_liquid.mat'])

end
