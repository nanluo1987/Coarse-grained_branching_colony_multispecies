%% Population dynamics
figure('position',[712   289   243   215])

% area(BiomassV)
nt = totalt / dt;
Frat = BiomassV ./ sum(BiomassV, 2); 
area((0 : nt) * dt, Frat)

color1 = [199,61,120]; color2 = [255,192,0]; color3 = [52,117,166];
colororder([color1; color2; color3]/255)
ylim([0 1]); xlim([0 totalt])
xlabel 'Time (h)'
ylabel 'Fraction'
set(gca, 'TickLength', [0,0])
set(gca, 'FontSize', 8)
set(gca, 'LabelFontSizeMultiplier', 10/8)
set(gca, 'xtick', 0:4:16)

% legend('WT','Cheater','Hyperswarmer')
% legend('location','northeastoutside')
% legend('box','off')

%% Multi days
foldername = [pwd '\results\'];
filename = dir([foldername '*_7d.mat']);
filename = {filename.name}';
nn = length(filename);

% figure('position',[455.6667  427.0000  178.6667   99.3333])
figure('position', [320.3333  345.6667  185.3333   82.6667])

for ii = 1 : nn

clf
load([foldername, filename{ii}])
nt = totalt / dt;
Frat = BiomassMultiDays ./ sum(BiomassMultiDays, 2); 
area(linspace(0, ndays, (nt + 1) * ndays + 1), Frat)

color1 = [199,61,120]; color2 = [255,192,0]; color3 = [52,117,166];
colororder([color1; color2; color3]/255)
ylim([0 1]); xlim([0 ndays])
xlabel 'Days'
ylabel 'Fraction'
set(gca, 'TickLength', [0,0])
set(gca, 'FontSize', 8)
set(gca, 'LabelFontSizeMultiplier', 10/8)
set(gca, 'xtick', 0 : ndays)
set(gca, 'ytick', [0 1])
set(gca, 'TitleFontWeight', 'bold')
set(gca, 'XGrid', 'on')
set(gca, 'Layer', 'top')
set(gca, 'GridAlpha', 0.5)

fig = gcf; fig.PaperPositionMode = 'auto';
print([foldername, filename{ii}(1 : end - 3) 'tif'],'-dtiff','-r800')

end