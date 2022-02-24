clear;
np = 1000;
load 'parameters_local_2_cell7_ID681_tuned2.mat'
Parameters_origin = Parameters;
perturbmag = 0.2; % magnitude of perturbation to parameters
ParametersMat = Parameters .* ones(1, np);
ParametersMat = ParametersMat .* (1 + perturbmag *(2 * rand(size(Parameters,1), np) - 1)); % perturb each parameter p within the [1-perturbmag, 1+perturbmag]*p range
ParametersMat = roundn(ParametersMat, -1);
ParametersMat(12, :) = roundn(rand(1, np) * 0.2 + 0.8, -2);
ParametersMat(13, :) = roundn(rand(1, np) * 0.2 + 0.8, -2);
ParametersMat(:, 1) = [Parameters_origin; 1; 0.9];
save('screen_new_parameters.mat', 'np', 'ParametersMat')

%% Round 1: WT/CT & WT/HS
clear; 
load 'screen_new_parameters.mat'
% parpool(8)

parfor ID = 1 : np
    
prefix = num2str(ID, '%04d');
Parameters = ParametersMat(:, ID);

for iter = 4 : 5
    
    initRatios = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1]; 
    initialRatio = initRatios(iter, :);   % initial ratio of all species
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species

    filename = ['results\screen_new\' prefix];
    Packed_for_screening(filename, initialFract, initialRatio, Parameters, iter)

end
disp(ID)

end

%% Round 2: WT/CT & WT/HS competition

% clear; 
load 'screen_new_parameters.mat'
initRatios = [1 0.01 0]; ra = 2; rb = 1; 
% initRatios = [1 0 0.01]; ra = 3; rb = 1; 
% parpool(8)
N0V = [1, 1.4, 1.6, 1.8]; 
% finalRatiosMat = zeros(length(N0V), np);
figure(2)
set(gcf,'position',[541.6667  375.0000  293.3333  232.0000])
load('results\screen_new\WC_varyN0.mat');
% finalRatioMat = ld.finalRatiosMat;

for ID = 901 : 1000
    
prefix = [num2str(ID, '%04d') '_WC'];

if exist(['results\screen_new\' num2str(ID, '%04d') '_04.jpg'],'file')
    
    finalRatioV = zeros(length(N0V), 1);
    
for iN = 1 : length(N0V)
    
    Parameters = ParametersMat(:, ID);
    Parameters(3)  = Parameters(3) * N0V(iN);
    Parameters(10) = Parameters(10) / N0V(iN);
    Parameters(11) = Parameters(11) / N0V(iN);

    initialRatio = initRatios;   % initial ratio of all species
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species

    filename = ['results\screen_new\' prefix];
    
    figure(1)
    SetParameters
    BranchingColonyMultispecies_Core
%     SaveFigure
    
    finalRatio = Biomass(ra) ./ Biomass(rb);
    finalRatiosMat(iN, ID) = finalRatio;
%     finalRatioV(iN) = finalRatio;
    
end

figure(2); clf; 
bar(N0V,finalRatiosMat(:, [1,ID])/0.01)
% bar(N0V, [finalRatiosMat(:, 1) finalRatioV] / 0.01)
set(gca, 'YScale', 'log'); ylim([1 inf])
saveas(gcf, [filename '.jpg'])
disp(ID)

end

save('results\screen_new\WC_varyN0.mat', 'N0V', 'ParametersMat', 'finalRatiosMat')

end

%% Round 3: iter 7

clear; 
load 'screen_new_parameters.mat'
% parpool(8)

for ID = 801 : 1000
    
prefix = num2str(ID, '%04d');
Parameters = ParametersMat(:, ID);

if exist(['results\screen_new\' num2str(ID, '%04d') '_WC.jpg'],'file')

for iter = 1
    
    initRatios = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1]; 
    initialRatio = initRatios(iter, :);   % initial ratio of all species
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species

    filename = ['results\screen_new\' prefix];
    Packed_for_screening(filename, initialFract, initialRatio, Parameters, iter)

end
disp(ID)

end

end

%% Round 4: competition

clear; 
load 'screen_new_parameters.mat'
% parpool(8)

for ID = 1 : 1000
    
prefix = num2str(ID, '%04d');
Parameters = ParametersMat(:, ID);

% if exist(['results\screen_new\' num2str(ID, '%04d') '_WC.jpg'],'file')
    
    finalratioV = zeros(4, 1);
    
for iter = 1 : 5
    
    CH_ratios = 10.^(linspace(-3,1,5))'; 
    initRatios = [zeros(length(CH_ratios),1), CH_ratios, ones(length(CH_ratios),1)]; 
    initialRatio = initRatios(iter, :);   % initial ratio of all species
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species

    filename = ['results\screen_new\' prefix];
    SetParameters
    figure(1)
    BranchingColonyMultispecies_Core
    
    finalRatio = Biomass(2) ./ Biomass(3);
    finalratioV(iter) = finalRatio;
    
end

figure(2); clf; 
plot(CH_ratios', finalratioV,'o-','linewidth', 2); hold on
plot(CH_ratios', [0.1848  0.3904   1.1171    0.3693    0.1325], 'o-','linewidth', 2)
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log'); ylim([0.01, 100])
saveas(gcf, [filename '_CH_compt.jpg'])
disp(ID)

% end

end

%% 
finished = zeros(1000, 1); 
for ID = 1 : 1000
    if exist(['results\screen_new\' num2str(ID, '%04d') '_04.jpg'],'file')
        finished(ID) = 1;
        if exist(['results\screen_new\' num2str(ID, '%04d') '_WH.jpg'],'file')
            finished(ID) = 2;
            if exist(['results\screen_new\' num2str(ID, '%04d') '_WC.jpg'],'file')
                finished(ID) = 3;
                if exist(['results\screen_new\' num2str(ID, '%04d') '_07.jpg'],'file')
                    finished(ID) = 4;
                    if exist(['results\screen_new\' num2str(ID, '%04d') '_01.jpg'],'file')
                        finished(ID) = 5;
                    end
                end
            end
        end
    end
end
figure; heatmap(reshape(finished, 50,20),'Colormap',hot)