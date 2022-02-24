%% One ratio
clear
initialRatio = [1, 1, 1];   % initial ratio of all species
initialFract = initialRatio / sum(initialRatio); % initial fraction of each species

SetParameters
BranchingColonyMultispecies_Core

%% Multiple ratios & varying N0
clear; 

initRatios = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1]; prefix = 'pyramid_744'; ra = 1; rb = 1;
% WC_ratios = 10.^(linspace(-3,1,41))'; initRatios = [ones(length(WC_ratios),1), WC_ratios, zeros(length(WC_ratios),1)]; prefix = 'WC'; ra = 2; rb = 1;
% CH_ratios = 10.^(linspace(-4, 4, 25))'; initRatios = [zeros(length(CH_ratios),1), CH_ratios, ones(length(CH_ratios),1)]; prefix = 'CH_C0x0.01'; ra = 2; rb = 3;
% WH_ratios = 10.^(linspace(-3,1,41))'; initRatios = [ones(length(WH_ratios),1), zeros(length(WH_ratios),1), WH_ratios]; prefix = 'WH'; ra = 3; rb = 1;
% initRatios = [1 1 1; 1 0.01 0.01; 0.5 0.5 0.01; 0.5 0.01 0.5]; prefix = '3sp'; ra = 1; rb = 1;
% initRatios = [1 1 1; 1 0.01 0.01]; prefix = '3sp'; ra = 1; rb = 1;

N0V = 1; 
% N0V = 1;
finalRatios_N0 = zeros(size(initRatios, 1), length(N0V));
% frontNutrient  = cell(3, 1); % #### recording nutrient

for iN = 1 : length(N0V)

fprintf('N0 = %f\n', N0V(iN))
load 'parameters_local_2_cell7_ID681_tuned2.mat'
% load 'parameters_local_2_cell6_ID744.mat'
Parameters(3)  = Parameters(3) * N0V(iN);
Parameters(10) = Parameters(10) / N0V(iN);
Parameters(11) = Parameters(11) / N0V(iN);

Output_Biomass = zeros(size(initRatios, 1), 3);
Output_Sizes   = zeros(size(initRatios, 1), 3);
Output_Biomass_Liq = zeros(size(initRatios, 1), 3);

for iter = 3
    
    fprintf('iter = %d\n', iter)
    initialRatio = initRatios(iter, :);   % initial ratio of all species
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species
    filename = ['results\' prefix '_N' num2str(iN)];
    SetParameters
    gs(3) = 1;
    gs = gs * 2;
    C0 = 0.8;
    save([filename '_parameters.mat'])
        
    filename = [filename '_' num2str(iter,'%02d')];
    BranchingColonyMultispecies_Core_new1
    SaveFigure
    save([filename '.mat'], 'BiomassV', 'C')
    % save([picname '_' num2str(iter,'%02d') '_' num2str(iN) '.mat'])
    Output_Biomass(iter, :) = Biomass;
%     Output_Sizes(iter, :) = Sizes;
    
%     LiquidCulture
%     Output_Biomass_Liq(iter, :) = Biomass;
%     save([filename '_liq.mat'], 'BiomassV', 'N', 'C')
    
    
end

finalRatio = Output_Biomass(:, ra) ./ Output_Biomass(:, rb);
finalRatios_N0(:, iN) = finalRatio;
save([filename '_results.mat'], 'initRatios', 'Output_Biomass', 'Output_Sizes', 'Output_Biomass_Liq')

end

save(['results\' prefix '_varyN0.mat'], 'N0V', 'finalRatios_N0')

%% Multiple cycles (days)

clear; 

ndays = 7;
initRatios = [1 1e-6 1e-6]; % ; 1 1e-4 1e-4; 1 0.01 0.01; 1 0.01 0; 1 0 0.01; 1 0.01 0; 1 0 0.01];
addRatios  = [0 0 0]; %       0 0 0;    0 0 0.01; 0 0.01 0];
prefix = 'evo_'; 
N0V = [0.8]; 
dayj = 4; % day # to add the third species

for iN = 1

fprintf('N0 = %f\n', N0V(iN))
load 'parameters_local_2_cell7_ID681_tuned2.mat'
Parameters(3)  = Parameters(3) * N0V(iN);
Parameters(10) = Parameters(10) / N0V(iN);
Parameters(11) = Parameters(11) / N0V(iN);

Output_Biomass = zeros(size(initRatios, 1), 3);
Output_Sizes   = zeros(size(initRatios, 1), 3);
Output_Biomass_Liq = zeros(size(initRatios, 1), 3);

for iter = 1 : size(initRatios, 1)
    
    initialRatio = initRatios(iter, :);   % initial ratio of all species
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species
    SetParameters
    if iN == 1 && iter == 1; save(['results\' prefix '_' num2str(ndays) 'd_parameters.mat']); end
    BiomassMultiDays = [];
    
    for iday = 1 : ndays
        BranchingColonyMultispecies_Core
        filename = ['results\' prefix '_N' num2str(N0V(iN)) '_' num2str(iter,'%02d') '_d' num2str(iday)];
        SaveFigure
        BiomassMultiDays = [BiomassMultiDays; BiomassV(1 : end - 1,:)]; 
        initialRatio = BiomassV(end, :);
        initialFract = initialRatio / sum(initialRatio);
        if iday == dayj
            initialFract = initialFract + addRatios(iter, :); 
            initialFract = initialFract / sum(initialFract);
        end
        fprintf('Day %d\n', iday)
    end
    BiomassMultiDays = [BiomassMultiDays; BiomassV(end, :)];
    
    save(['results\' prefix '_N' num2str(N0V(iN)) '_' num2str(iter,'%02d') '_' num2str(ndays) 'd.mat'], 'BiomassMultiDays','totalt','dt','ndays')
    Output_Biomass(iter, :) = Biomass;
    Output_Sizes(iter, :) = Sizes;
    fprintf('iter = %d\n', iter)
    
end

end

%% Test: HS lower growth rate
clear; 

pslow = 1;
% WH_ratios = 0.01; initRatios = [ones(length(WH_ratios),1), zeros(length(WH_ratios),1), WH_ratios]; prefix = 'WH'; ra = 3; rb = 1;
WC_ratios = 0.01; initRatios = [ones(length(WC_ratios),1), WC_ratios, zeros(length(WC_ratios),1)]; prefix = 'WC'; ra = 2; rb = 1;
N0V = 1.7; finalRatios_N0 = zeros(size(initRatios, 1), length(N0V));

for ipslow = 1 : length(pslow)
    
for iN = 1 : length(N0V)

fprintf('N0 = %f\n', N0V(iN))
load 'parameters_local_2_cell7_ID681_tuned2.mat'
Parameters(10) = 0.6;
Parameters(3)  = Parameters(3) * N0V(iN);
Parameters(10) = Parameters(10) / N0V(iN);
Parameters(11) = Parameters(11) / N0V(iN);

Output_Biomass = zeros(size(initRatios, 1), 3);
Output_Sizes   = zeros(size(initRatios, 1), 3);
Output_Biomass_Liq = zeros(size(initRatios, 1), 3);

for iter = 1 : size(initRatios, 1)
    
    initialRatio = initRatios(iter, :);   % initial ratio of all species
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species
    filename = ['Lower HS growth\' prefix '_N' num2str(iN)];
    SetParameters
    aCs_act(3) = pslow(ipslow) * Parameters(4); % ### lower HS growth rate ###
%     gs(3) = 2.5 * Parameters(5);   
    save([filename '_parameters.mat'])
    
    filename = [filename '_' num2str(iter,'%02d')]; % '_' num2str(pslow(ipslow),'%.2f')];
    BranchingColonyMultispecies_Core
    SaveFigure
    save([filename '.mat'], 'BiomassV', 'C')
    % save([picname '_' num2str(iter,'%02d') '_' num2str(iN) '.mat'])
    Output_Biomass(iter, :) = Biomass;
    Output_Sizes(iter, :) = Sizes;
    
    LiquidCulture
    Output_Biomass_Liq(iter, :) = Biomass;
    save([filename '_liq.mat'], 'BiomassV', 'N', 'C')
    fprintf('iter = %d\n', iter)
    
end

finalRatio = Output_Biomass(:, ra) ./ Output_Biomass(:, rb);
finalRatios_N0(:, iN) = finalRatio;
save([filename '_results.mat'], 'initRatios', 'Output_Biomass', 'Output_Sizes', 'Output_Biomass_Liq')

end

save(['Lower HS growth\' prefix '_varyN0_upper0.6.mat'], 'N0V', 'finalRatios_N0')

end