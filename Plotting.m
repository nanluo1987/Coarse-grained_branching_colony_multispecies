if mod(i, 10) == 0
    ind = 1 : 2 : nx;
    
    figure(1)
    % Plot each species
    subplot(2, 3, j)
        hold off; pcolor(xx(ind, ind), yy(ind, ind), C{j}(ind, ind));
        shading interp; axis equal;
        axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
        colorbar
        set(gca,'YTick',[], 'XTick',[])
        plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
        plot(Tipx{j}(aCs_tip < max(aCs_act), ib), Tipy{j}(aCs_tip < max(aCs_act), ib), '.r', 'markersize', 5) % plot swarming tips
        title(speciesName{j})
        drawnow

    % Plot all species
    if j == find(initialRatio,1,'last')

    Ctotal = 5*C{1} + C{2} + C{3};
    p1 = 5*C{1}./Ctotal; p1(isnan(p1)) = 0; 
    p2 = C{2}./Ctotal; p2(isnan(p2)) = 0;
    p3 = C{3}./Ctotal; p3(isnan(p3)) = 0;
    
    ind = 1 : 2 : nx;
%     color1 = [199,61,120]; color2 = [255,192,0]; color3 = [52,117,166];
    color1 = [251,69,199]; color2 = [255,255,0]; color3 = [17,207,226];
    subplot(2, 3, 4) % total cell density
%         hold off; pcolor(xx(ind, ind), yy(ind, ind), Ctotal(ind, ind));
%         shading interp; axis equal;
%         axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
%         colorbar
%         set(gca,'YTick',[], 'XTick',[])
%         plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
%         title(['Time = ' num2str(i * dt)])
    subplot(2, 3, 5) % show each species by color
        ColorMap = MarkMixing_3color(color1, color2, color3, p1, p2, p3);
        ColorMap(isnan(ColorMap)) = 0;
        hold off; surf(xx(ind, ind), yy(ind, ind), ones(size(xx(ind, ind))), ColorMap(ind, ind, :))
        view([0, 0, 1]); shading interp; axis equal; box on
        axis([-L/2 L/2 -L/2 L/2]);
        set(gca,'YTick',[], 'XTick',[])
        title(['Time = ' num2str(i * dt)])
    subplot(2, 3, 6) % line graph of cell densities
        yyaxis left; hold off
        mid = (nx + 1) / 2;
        plot(x(mid:end), C{1}(mid:end,mid), '-', 'color', color1/255, 'linewidth', 2); hold on
        plot(x(mid:end), C{2}(mid:end,mid), '-', 'color', color2/255, 'linewidth', 2);
        plot(x(mid:end), C{3}(mid:end,mid), '-', 'color', color3/255, 'linewidth', 2);
        plot(x(mid:end), Ctotal(mid:end,mid), 'k-', 'linewidth', 2)
        ylabel 'Cell density';
        yyaxis right; hold off
        plot(x(mid:end), N(mid:end,mid), '-', 'color', [0.7,0.7,0.7], 'linewidth', 2); ylim([0 N0])
        xlabel 'Distance from center'
%         hold off; pcolor(xx(ind, ind), yy(ind, ind), hs{3}(ind, ind));
%         shading interp; axis equal;
%         axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
%         colorbar
%         set(gca,'YTick',[], 'XTick',[])
%         title('Nutrient')
    drawnow
    end
    
%     figure(8)
%     csize = mean(sqrt(Tipx{j}(1:nn,ib).^2+Tipy{j}(1:nn,ib).^2)); % colony size
%     hstip = interp2(xx, yy, hs{j}, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib)); 
%     swarmingtips = sum(hstip > 0); % number of swarming tip
%     gs_fract = mean(gs(j) ./ (gs(j) * 1 + aCs_tip .* sum(dE, 2))); % swimming motility fraction
%     if j == 1; cc = 'or';
%     elseif j== 3; cc = 'ob'; end
%     subplot 221; plot(i, mean(dl), cc); hold on; 
%     subplot 222; plot(i, csize, cc); hold on; 
%     subplot 223; plot(i, swarmingtips, cc); hold on; 
%     subplot 224; plot(i, gs_fract, cc); hold on; 
%     drawnow
%     if j == find(initialRatio,1,'first'); idxN = idxN + 1; end
%     tipN = mean(interp2(xx, yy, N, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib)));   
%     frontNutrient{j}(iN, idxN) = tipN;
    
end