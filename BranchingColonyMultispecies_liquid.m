clear
clf

for iter = 1 : 5
    
figure(1)

%% Parameters

totalt = 14;    % total time
dt     = 0.02;  % time step

speciesName = {'WT','Cheater','Hyperswarmer'}; % name of each species
% other vectors will follow the same order
% pyramid_initRatios = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; 0.9, 0.1, 0; 0.99 0.01 0];
% pyramid_initRatios = [0.99 0.01 0; 0.9 0.1 0; 1 1 0];
pyramid_initRatios = [1 0 0; 100 1 0; 10 1 0; 1 1 0; 1 10 0];
initialRatio = pyramid_initRatios(iter, :);   % initial ratio of all species
initialFract = initialRatio / sum(initialRatio); % initial fraction of each species

bN = 25;         % nutrient consumption rate
DN = 7;          % nutrient diffusivity
N0 = 20;         % initial nutrient conc.
 
aCs_act = [1, 1.2, 1] * 0.75;      % cell growth rate of each species
gs  = [1, 1, 2] * 1;            % swimming motility
hs_act  = [1, 0, 0.9] * 8;        % swarming motility coefficients
mu = 0.4;

N_upper = 15; % upper bound of nutrient for swarming
N_lower = 6; % lower bound of nutrient for swarming

%% Initialization
C0  = 8;   % initial cell density
nt = totalt / dt;  % number of time steps

C  = C0 /100 * initialFract;
N  = N0;

Cv = zeros(nt + 1, 3); Cv(1, :) = C; % Cell density
Nv = zeros(nt + 1, 1); N(1, :)  = N; % Nutrient

aCs = zeros(3, 1);

for i = 1 : nt
    
    fN = N / (N + 1) * (1 - sum(C));
    if N > N_upper || N < N_lower
        aCs = max(aCs_act) * ones(1, 3);
    else
        aCs = aCs_act;
    end
    
    C = C + aCs .* C * fN * dt;     
    N = N - bN * fN .* sum(aCs .* C) * dt;
    
    Cv(i + 1, :) = C;
    Nv(i + 1, :) = N;
            
end

subplot 121; 
yyaxis left; cla; hold on; 
plot(0 : dt : totalt, Cv, 'linewidth', 2)
ylabel 'Cell density';
yyaxis right; cla; hold on
plot(0 : dt : totalt, Nv, '-', 'color', [0.7,0.7,0.7], 'linewidth', 2); ylim([0 N0])
ylabel 'Nutrient'
xlabel 'Time'
drawnow

subplot 122
hold on
plot(initialRatio(2)/initialRatio(1), C(2)/C(1), 'o')
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

end
