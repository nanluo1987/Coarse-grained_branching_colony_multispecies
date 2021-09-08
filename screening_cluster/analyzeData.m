clf; clear
load('output.mat')
pyramid_initRatios = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1];
     
DesiredSizes = [40,0,0,0,5,0,0,0,45,27,27,0,25,0,40,0,18,27,10,18,27]' * ones(1, 10000);
inonzero = DesiredSizes(:,1) > 0;
SizeScore = abs(DesiredSizes(inonzero,:) - OutputMat_Sizes(inonzero,:)) ./ DesiredSizes(inonzero,:);
SizeScore_rel = SizeScore - SizeScore(:, 1);

Ratios = zeros(5, 10000);
Ratios(1, :) = OutputMat_Sizes(10, :) ./ OutputMat_Sizes(11, :); % WT:CT in WT/CT
Ratios(2, :) = OutputMat_Sizes(13, :) ./ OutputMat_Sizes(15, :); % WT:HS in WT/HS
Ratios(3, :) = OutputMat_Sizes(17, :) ./ OutputMat_Sizes(18, :); % CT:HS in CT/HS
Ratios(4, :) = OutputMat_Sizes(19, :) ./ OutputMat_Sizes(20, :); % WT:CT in WT/CT/HS
Ratios(5, :) = OutputMat_Sizes(20, :) ./ OutputMat_Sizes(21, :); % CT:HS in WT/CT/HS
DesiredRatios = [DesiredSizes(10,1)./DesiredSizes(11,1);
                 DesiredSizes(13,1)./DesiredSizes(15,1);
                 DesiredSizes(17,1)./DesiredSizes(18,1);
                 DesiredSizes(19,1)./DesiredSizes(20,1);
                 DesiredSizes(20,1)./DesiredSizes(21,1)] * ones(1, 10000);
RatioScore = abs(DesiredRatios - Ratios) ./ DesiredRatios;
RatioScore_rel = RatioScore - RatioScore(:, 1);

i_bottomline = Ratios(1, :) > 0.9 & Ratios(2, :) < 1 & ...
               Ratios(3, :) < 1 & Ratios(4, :) < 1 & Ratios(5, :) < 1 ...
               & OutputMat_Biomass(1, :) > sum(OutputMat_Biomass(end-2:end,:),1);
keyScore = [SizeScore([1,3],:); RatioScore];
% i_30perc = sum(keyScore < 0.3, 1) >= 6;
i_selected = i_bottomline ...
           & keyScore(1,:) < 0.35 ... % WT size
           & keyScore(2,:) < 0.3 ... % HS size
           & keyScore(3,:) < 0.05 ... % WT/CT ratio in WT-CT
           & keyScore(6,:) < 0.38; ... % WT/CT ratio in 3-sp

heatmap(find(i_selected),{'WT size','HS size','WT/CT ratio in WT-CT','WT/HS ratio in WT-HS','CT/HS ratio in CT-HS','WT/CT ratio','CT/HS ratio'},...
        keyScore(:, i_selected)); colormap 'parula'
IDV = find(i_selected); %IDV = [1, IDV];

for ID = IDV
    disp(ID)
    BranchingColonyMultispecies_reprod
end
    

% WT size
% iter = 1; 
% plot(OutputMat_Sizes((iter - 1) * 3 + 1, selected), 'o'); hold on
% 
% % HS size
% iter = 3; 
% plot(OutputMat_Sizes((iter - 1) * 3 + 3, selected), 'o')
% 
% % CT size
% iter = 2; 
% plot(OutputMat_Sizes((iter - 1) * 3 + 2, selected), 'o')