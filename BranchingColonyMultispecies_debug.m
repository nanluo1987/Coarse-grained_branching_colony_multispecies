load('189.mat')

for i = 190 : nt

    % -------------------------------------
    % Nutrient distribution

    fN = N ./ (N + 1) .* (1 - (C{1} + C{2} + C{3}));
    dN = - bN * fN .* (aCs(1) * C{1} + aCs(2) * C{2} + aCs(3) * C{3}); % Nutrient consumption
    N  = N + dN * dt;
    NV = MatV1N \ (N * MatU1N); N = (MatV2N * NV) / MatU2N; % Nutrient diffusion

    for j = 1 : 3  % j: index of species

    % -------------------------------------
    % Cell growth

    C{j} = C{j} + aCs(j) * fN .* C{j} * dt;

    % -------------------------------------
    % Cell allocation and branch extension

    if mod(i, dt_updatebranch / dt) == 0 && initialFract(j) > 0
        
        ib = i / (dt_updatebranch / dt) + 2;
        Tipx{j}(:, ib) = Tipx{j}(:, max(1, ib - 1));
        Tipy{j}(:, ib) = Tipy{j}(:, max(1, ib - 1));
        dBiomass = (C{j} - C_pre{j}) * dx * dy; % total cell growth

        % compute the amount of biomass accumulation in each branch
        BranchDomainSum = cat(3, BranchDomain{:,j});
        BranchDomainSum = sum(BranchDomainSum, 3);
        nn = ntips;
        for k = 1 : nn
            branchfract = 1 ./ (BranchDomainSum .* BranchDomain{k,j});
            branchfract(isinf(branchfract)) = 0;
            dE(k,j) = sum(sum(dBiomass .* sparse(branchfract)));
        end
        
        % update width and density with time. 
        currFracts = [sum(C{1}, 'all') sum(C{2}, 'all') sum(C{3}, 'all')];
        currFracts = currFracts ./ sum(currFracts); %single vector, the below is spatial grid. 
        Density = currFracts * Densities';
        Width   = currFracts * Widths';
        
        % gamma (expansion efficiency) = swimming efficiency + swarming efficiency
        % calculate the gamma at each branch tip from the local composition
        tipC = [interp2(xx, yy, C{1}, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib)) ...
                interp2(xx, yy, C{2}, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib)) ...
                interp2(xx, yy, C{3}, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib))];
        tipFract = tipC ./ sum(tipC, 2);
        gammas = gs(j) + sum(hs .* tipFract, 2);
        
        % extension rate of each branch  
        dl = gammas .* dE(1:nn,j) ./ Width;
        if i == 0; dl = 0.5; end

        [Tipx{j}, Tipy{j}] = tiptracking(Tipx, Tipy, ib, dl, theta, delta, nn, xx, yy, N, j, noiseamp);
        
        % Bifurcation
        R = 3/2 / Density;  % a branch will bifurcate if there is no other branch tips within the radius of R
        TipxNew = Tipx{j}; TipyNew = Tipy{j};
        thetaNew = theta(:,j); dlNew = dl;
        BranchDomainNew = BranchDomain(:,j);
        for k = 1 : nn
            dist2othertips = sqrt((TipxNew(:,ib) - Tipx{j}(k,ib)) .^ 2 + (TipyNew(:,ib) - Tipy{j}(k,ib)) .^ 2);
            dist2othertips = sort(dist2othertips);
            if dist2othertips(2) > R
                nn = nn + 1;
                TipxNew(nn,:) = Tipx{j}(k,:);
                TipyNew(nn,:) = Tipy{j}(k,:);
                for jk = 1 : 3
                    if jk ~= j && size(Tipx{jk}, 1) < nn
                        Tipx{jk}(nn,:) = Tipx{jk}(k,:);
                        Tipy{jk}(nn,:) = Tipy{jk}(k,:);
                        BranchDomain{nn,jk} = BranchDomain{k,jk};
                    end
                end
                TipxNew(nn,ib) = Tipx{j}(k,ib) + dl(k) * sin(theta(k,j) + 0.5 * pi); % splitting the old tip to two new tips
                TipyNew(nn,ib) = Tipy{j}(k,ib) + dl(k) * cos(theta(k,j) + 0.5 * pi);
                TipxNew(k,ib) = TipxNew(k,ib) + dl(k) * sin(theta(k,j) - 0.5 * pi);
                TipyNew(k,ib) = TipyNew(k,ib) + dl(k) * cos(theta(k,j) - 0.5 * pi);
                dlNew(nn) = dl(k) / 2;
                dlNew(k)  = dl(k) / 2;
                thetaNew(nn) = theta(k,j);
                BranchDomainNew{nn} = BranchDomain{k,j};
            end
        end
        ntips = nn;
        if nn > size(Tipx{j}, 1)
            nz = nn - size(Tipx{j}, 1);
            theta = [theta; zeros(nz, 3)];
        end
        Tipx{j} = TipxNew; Tipy{j} = TipyNew;
        theta(:,j) = thetaNew; dl = dlNew;
        BranchDomain(:,j) = BranchDomainNew;    
        
        

        % Fill the width of the branches
        for k = 1 : nn
            d = sqrt((Tipx{j}(k,ib) - xx) .^ 2 + (Tipy{j}(k,ib) - yy) .^ 2);
            P{j}(d <= Width/2) = 1;
            BranchDomain{k,j} = BranchDomain{k,j} | (d <= Width/2);
        end

        % relocate dBiomass
        Capacity = 1 - C_pre{1} - C_pre{2} - C_pre{3};    % remaining cell capacity
        Capacity(P{j} == 0) = 0;   % no capacity outside the colony
        frac_relo = 1;
        C_relo = frac_relo * sum(dBiomass(:)) / sum(Capacity(:)) * Capacity / (dx * dy);
        C{j} = C_pre{j} + C_relo + (1 - frac_relo) * dBiomass / (dx * dy);
        C_pre{j} = C{j};

        % Plot each species
        ind = 1 : 2 : nx;
        subplot(2, 3, j)
            hold off; pcolor(xx(ind, ind), yy(ind, ind), C{j}(ind, ind));
            shading interp; axis equal;
            axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
            colorbar
            set(gca,'YTick',[], 'XTick',[])
            plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
            title(speciesName{j})
            drawnow
        
        if j == 2
        % Plot all species
        Ctotal = C{1} + C{2} + C{3};
        p1 = C{1}./Ctotal; p1(isnan(p1)) = 0;
        p2 = C{2}./Ctotal; p2(isnan(p2)) = 0;
        p3 = C{3}./Ctotal; p3(isnan(p3)) = 0;
        ind = 1 : 2 : nx;
        color1 = [199,61,120]; color2 = [255,192,0]; color3 = [52,117,166];
        subplot(2, 3, 4) % total cell density
            hold off; pcolor(xx(ind, ind), yy(ind, ind), Ctotal(ind, ind));
            shading interp; axis equal;
            axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
            colorbar
            set(gca,'YTick',[], 'XTick',[])
            plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
            title(['Time = ' num2str(i * dt)])
        subplot(2, 3, 5) % show each species by color
            ColorMap = MarkMixing_3color(color1, color2, color3, p1, p2, p3);
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
        drawnow
        end
        
    end

    end
    
    % Growth stops when approaching edges
%     TipR = sqrt(Tipx.^2 + Tipy.^2);
%     if max(TipR(:)) > 0.9 * L/2
%         break
%     end

end