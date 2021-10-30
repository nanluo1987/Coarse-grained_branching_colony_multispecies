clf; clear
% load('output.mat'); njobs = 20000;
load('output_local_2.mat'); njobs = 2000;
iP = 7; 
ParametersMat = ParametersCell{iP};
OutputMat_Sizes = OutputCell_Sizes{iP};
OutputMat_Biomass = OutputCell_Biomass{iP};

pyramid_initRatios = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; 0.9 0.1 0; 0.99 0.01 0; 0.1 0.9 0];
% output formats:
% OutputMat = (3 * n_initRatios) * njobs

%% Primary criteria

% Biomass of 3-species colony < biomass of WT
idx1 = sum(OutputMat_Biomass(19:21,:),1) < OutputMat_Biomass(1,:);
% Radius of WT < radius of HS in WT+HS colonies
idx2 = OutputMat_Sizes(13,:) < OutputMat_Sizes(15,:);
% Radius of CT < radius of HS in CT+HS colonies
idx3 = OutputMat_Sizes(17,:) < OutputMat_Sizes(18,:);
% Radius of WT < radius of CT < radius of HS in 3-species colonies
idx4 = OutputMat_Sizes(19,:) < OutputMat_Sizes(20,:) & OutputMat_Sizes(20,:) < OutputMat_Sizes(21,:);
% CT outcompetes WT at starting proportion of 0.1 or 0.01
idx5 = OutputMat_Biomass(23,:) ./ (OutputMat_Biomass(22,:) + OutputMat_Biomass(23,:)) > 0.1 ...
     & OutputMat_Biomass(26,:) ./ (OutputMat_Biomass(26,:) + OutputMat_Biomass(25,:)) > 0.01;
 
selected_pri = idx1 & idx2 & idx3 & idx4 & idx5;
disp(sum(selected_pri))

finalFrat_8 = OutputMat_Biomass(23,:) ./ (OutputMat_Biomass(22,:) + OutputMat_Biomass(23,:));
finalFrat_9 = OutputMat_Biomass(26,:) ./ (OutputMat_Biomass(26,:) + OutputMat_Biomass(25,:));

%% Secondary criteria

DesiredSizes = [40,0,0,0,5,0,0,0,45,27,27,0,25,0,40,0,18,27,10,18,27]' * ones(1, njobs);
inonzero = DesiredSizes(:,1) > 0;
SizeScore = abs(DesiredSizes(inonzero,:) - OutputMat_Sizes(inonzero,:)) ./ DesiredSizes(inonzero,:);
% SizeScore_rel = SizeScore - SizeScore(:, 1);

Ratios = zeros(5, njobs);
Ratios(1, :) = OutputMat_Sizes(10, :) ./ OutputMat_Sizes(11, :); % WT:CT in WT/CT
Ratios(2, :) = OutputMat_Sizes(13, :) ./ OutputMat_Sizes(15, :); % WT:HS in WT/HS
Ratios(3, :) = OutputMat_Sizes(17, :) ./ OutputMat_Sizes(18, :); % CT:HS in CT/HS
Ratios(4, :) = OutputMat_Sizes(19, :) ./ OutputMat_Sizes(20, :); % WT:CT in WT/CT/HS
Ratios(5, :) = OutputMat_Sizes(20, :) ./ OutputMat_Sizes(21, :); % CT:HS in WT/CT/HS
DesiredRatios = [DesiredSizes(10,1)./DesiredSizes(11,1);
                 DesiredSizes(13,1)./DesiredSizes(15,1);
                 DesiredSizes(17,1)./DesiredSizes(18,1);
                 DesiredSizes(19,1)./DesiredSizes(20,1);
                 DesiredSizes(20,1)./DesiredSizes(21,1)] * ones(1, njobs);
RatioScore = abs(DesiredRatios - Ratios) ./ DesiredRatios;
% RatioScore_rel = RatioScore - RatioScore(:, 1);

% selected_sec = SizeScore(1,:) < 0.35 ... % WT size
%              & SizeScore(2,:) < 2 ...  % CT size
%              & SizeScore(3,:) < 0.3 ... % HS size 
%              & RatioScore(1,:) < 0.05 ... % WT/CT ratio in WT-CT
%              & RatioScore(4,:) < 0.38; ... % WT/CT ratio in 3-sp
selected_sec = selected_pri ...
             & OutputMat_Sizes(1,:) < 50 & OutputMat_Sizes(1,:) > 25 ... % WT size
             & OutputMat_Sizes(5,:) < 15 & OutputMat_Sizes(5,:) > 0 ... % CT size
             & OutputMat_Sizes(9,:) < 65 & OutputMat_Sizes(9,:) > 30 ... % HS size
             & Ratios(1,:) < 1.2 ... % WT/CT ratio in WT-CT
             & Ratios(4,:) < 0.9 ... % WT/CT ratio in 3-sp
             & Ratios(5,:) < 0.9 ... % CT/HS ratio in 3-sp
             & finalFrat_8 > 0.3 ... % Final fration of [0.9 0.1 0]
             & finalFrat_9 > 0.1; ... % Final fration of [0.99 0.01 0]
             
IDV = find(selected_sec); %IDV = [1, IDV];
disp(IDV); clf; 
subplot 121
heatmap(find(selected_sec),{'WT size','CT size','HS size', 'WT/CT ratio in WT-CT','WT/HS ratio in WT-HS','CT/HS ratio in CT-HS','WT/CT ratio','CT/HS ratio'},...
        [SizeScore(1:3, selected_sec); RatioScore(:, selected_sec)]); colormap 'parula'
subplot 122
plot(finalFrat_8(selected_sec), 'o'); hold on
plot(finalFrat_9(selected_sec), 'o')

%%
selected_sec = selected_pri ...
             & OutputMat_Sizes(1,:) < 52 & OutputMat_Sizes(1,:) > 38 ... % WT size
             & OutputMat_Sizes(5,:) < 12 & OutputMat_Sizes(5,:) > 0 ... % CT size
             & OutputMat_Sizes(9,:) < 52 & OutputMat_Sizes(9,:) > 38 ... % HS size
             & OutputMat_Sizes(11,:) > 15 ... % CT size in WT-CT
             & Ratios(1,:) > 0.88 ... % WT/CT ratio in WT-CT
             & Ratios(4,:) < 0.8 ... % WT/CT ratio in 3-sp
             & Ratios(5,:) < 0.85 & Ratios(5,:) > 0.6 ... % CT/HS ratio in 3-sp
             & finalFrat_8 > 0.5 ... % Final fration of [0.9 0.1 0]
             & finalFrat_9 > 0.2 & finalFrat_9 < 0.4; ... % Final fration of [0.99 0.01 0]

IDV = find(selected_sec); %IDV = [1, IDV];
disp(IDV); clf; 
disp(length(IDV))

subplot 911; bar(OutputMat_Sizes(1,[1,IDV])); ylabel 'WT'
subplot 912; bar(OutputMat_Sizes(5,[1,IDV])); ylabel 'CT'
subplot 913; bar(OutputMat_Sizes(9,[1,IDV])); ylabel 'HS'
subplot 914; bar(Ratios(1,[1,IDV])); ylabel 'WT/CT' % WT/CT ratio in WT-CT
subplot 915; bar(Ratios(4,[1,IDV])); ylabel 'WT/CT' % WT/CT ratio in 3-sp
subplot 916; bar(Ratios(5,[1,IDV])); ylabel 'CT/HS'
subplot 917; bar(finalFrat_8(:,[1,IDV])); ylabel 'Frat 0.1'
subplot 918; bar(finalFrat_9(:,[1,IDV])); ylabel 'Frat 0.01'
subplot 919; bar(OutputMat_Sizes(11,[1,IDV])); ylabel 'WT-CT size'


%%
IDV = 73;
for ID = IDV(1:end)
    disp(ID)
    picname = ['results\Local screening 2\refine\7-' num2str(ID,'%05d')];
%     picname = ['results\Local screening 2\' num2str(IPs(ID)) '-' num2str(IDs(ID),'%05d')];
%     if ~exist([picname '_01.jpg'], 'file')
    BranchingColonyMultispecies_reprod
%     end
end

%%

bestID = [00982 1485 769 191 991 141 1314];
ParameterSelected = zeros(11, 7);
OutputMat_Sizes_Selected = zeros(30, 7);
OutputMat_Biomass_Selected = zeros(30, 7);
for iP = 1 : 7  
    ParametersMat = ParametersCell{iP};
    OutputMat_Sizes = OutputCell_Sizes{iP};
    OutputMat_Biomass = OutputCell_Biomass{iP};
    ParameterSelected(:, iP) = ParametersMat(:, bestID(iP));
    OutputMat_Sizes_Selected(:, iP) = OutputMat_Sizes(:, bestID(iP));
    OutputMat_Biomass_Selected(:, iP) = OutputMat_Biomass(:, bestID(iP));
end
save('selectedParameters_local1.mat','ParameterSelected','OutputMat_Sizes_Selected','OutputMat_Biomass_Selected')

%% merge cells
load('output_local_2.mat'); njobs = 2000;
ParametersMat = [];
OutputMat_Sizes = [];
OutputMat_Biomass = [];
IPs = [];
IDs = [];

for iP = 1 : 7
    ParametersMat = [ParametersMat ParametersCell{iP}];
    OutputMat_Sizes = [OutputMat_Sizes OutputCell_Sizes{iP}];
    OutputMat_Biomass = [OutputMat_Biomass OutputCell_Biomass{iP}];
    IPs = [IPs; iP * ones(njobs,1)];
    IDs = [IDs; (1 : njobs)'];
end
njobs = length(IDs);

%% Get IDs
foldername = [pwd '\results\Local screening 2\'];
filename = dir([foldername '*_01.jpg']);
filename = {filename.name};
nn = length(filename);

IDV = zeros(1, nn);
for i = 1 : nn
    iP = str2double(filename{i}(1));
    ID = str2double(filename{i}(3:7));
    IDV(i) = find(IDs==ID & IPs == iP);
end

disp(IDV); clf; 
disp(length(IDV))

subplot 911; bar(OutputMat_Sizes(1,IDV)); ylabel 'WT'
subplot 912; bar(OutputMat_Sizes(5,IDV)); ylabel 'CT'
subplot 913; bar(OutputMat_Sizes(9,IDV)); ylabel 'HS'
subplot 914; bar(Ratios(1,IDV)); ylabel 'WT/CT' % WT/CT ratio in WT-CT
subplot 915; bar(Ratios(4,IDV)); ylabel 'WT/CT' % WT/CT ratio in 3-sp
subplot 916; bar(Ratios(5,IDV)); ylabel 'CT/HS'
subplot 917; bar(finalFrat_8(:,IDV)); ylabel 'Frat 0.1'
subplot 918; bar(finalFrat_9(:,IDV)); ylabel 'Frat 0.01'
subplot 919; bar(OutputMat_Sizes(11,IDV)); ylabel 'WT-CT size'

xticks(1:nn)
xticklabels(IPs(IDV))