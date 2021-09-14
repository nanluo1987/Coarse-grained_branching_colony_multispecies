figure(3); clf
colors = [199,61,120; 255,192,0; 52,117,166];

for ii = 1 : 14
    
    for j = 1 : 3
        
    load(['D:\Temp\111_' num2str(ii) '{' num2str(j) '}.mat'])
    TipR = sqrt(Tipx{j}(:,ib).^2 + Tipy{j}(:,ib).^2);
    
    subplot 131
    plot(ii, gammas * 0.06, 'o', 'color', colors(j, :)/255); hold on
    subplot 131
    plot(ii, gammas * 0.06+dE(1:nn,j) ./ Width, 'o', 'color', colors(j, :)/255); hold on
    subplot 132
    plot(ii, dE(1:nn,j) ./ Width, 'o', 'color', colors(j, :)/255); hold on
    subplot 133
    plot(ii, TipR, 'o', 'color', colors(j, :)/255); hold on
    drawnow
    
    end
    
end