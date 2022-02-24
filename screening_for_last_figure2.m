%% Multiple ratios & varying N0
clear; 
np = 6 * 6 * 6 + 1;
% WC_ratios = [0.01,1]'; initRatios = [ones(length(WC_ratios),1), WC_ratios, zeros(length(WC_ratios),1)]; prefix = 'WC'; ra = 2; rb = 1;
WH_ratios = 0.01; initRatios = [ones(length(WH_ratios),1), zeros(length(WH_ratios),1), WH_ratios]; prefix = 'screening_WH'; ra = 3; rb = 1;
% initRatios = [1 1 1; 1 0.01 0.01; 0.5 0.5 0.01; 0.5 0.01 0.5]; prefix = '3sp'; ra = 1; rb = 1;
% initRatios = [1 1 1; 1 0.01 0.01]; prefix = '3sp'; ra = 1; rb = 1;
N0V = [1.68]; 
finalRatiosMat = zeros(length(N0V), np);
NewParameters = zeros(13, np);

[aa,bb,cc] = meshgrid(linspace(1.5,3,6),linspace(0.8,1,6),linspace(0.8,1,6));
VarParameters = [aa(:) bb(:) cc(:)];
VarParameters = [2.3, 1, 0.9; VarParameters];

for ip = 1 : np
    
    for iN = 1 : length(N0V)

    fprintf('N0 = %f\n', N0V(iN))
    load 'parameters_local_2_cell7_ID681_tuned2.mat'
    Parameters(8)  = VarParameters(ip, 1);
    Parameters(12) = VarParameters(ip, 2);
    Parameters(13) = VarParameters(ip, 3);
    
    Parameters(3)  = Parameters(3) * N0V(iN);
    Parameters(10) = Parameters(10) / N0V(iN);
    Parameters(11) = Parameters(11) / N0V(iN);

    Output_Biomass = zeros(size(initRatios, 1), 3);
    Output_Sizes   = zeros(size(initRatios, 1), 3);
    Output_Biomass_Liq = zeros(size(initRatios, 1), 3);

    initialRatio = initRatios;   % initial ratio of all species
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species
%     filename = ['results\' prefix '_N' num2str(iN)];
    SetParameters
    NewParameters(:, ip) = Parameters;

%     filename = [filename '_' num2str(iter,'%02d')];
    BranchingColonyMultispecies_Core
%     SaveFigure    
    
    finalRatio = Biomass(ra) ./ Biomass(rb);
    finalRatiosMat(iN, ip) = finalRatio;

    end

fprintf('ip = %d\n', ip)
save(['results\' prefix '_3.mat'], 'N0V', 'NewParameters', 'finalRatiosMat')

end

