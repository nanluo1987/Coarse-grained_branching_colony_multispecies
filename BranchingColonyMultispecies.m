clear

%% Parameters
L      = 90;    % domain size
totalt = 16;    % total time
dt     = 0.02;  % time step
nx     = 1001; ny = nx; % number of nodes

speciesName = {'WT','Cheater','Hyperswarmer'}; % name of each species
% other vectors will follow the same order

initialRatio = [1, 0, 0];   % initial ratio of all species
initialFract = initialRatio / sum(initialRatio); % initial fraction of each species

bN = 25;             % nutrient consumption rate
DN = 7;              % nutrient diffusivity
N0 = 8/1.2;         % initial nutrient conc.

aCs = [1, 1.2, 1] * 1.5/2;      % cell growth rate of each species
gs = 0.4*[1 1 2];               % swimming motility
hs = [1, 0, 0.9] * 4.4;        % swarming motility coefficients

% branch density & width of single-species colonies
Densities = [0.14, 0.14, 0.2];
Widths    = [5, 5, 20];

% branch density & width of mixed colonies
Density = initialFract * Densities';
Width   = initialFract * Widths';

r0  = 5;     % initial radius
C0  = 8;   % initial cell density

noiseamp = 0;        % noise amplitude of branch direction
dt_updatebranch = 10 * dt;  % time step for updating branch locations

%% Initialization

nt = totalt / dt;  % number of time steps
dx = L / (nx - 1); dy = dx;
x  = linspace(-L/2, L/2, nx);
y  = linspace(-L/2, L/2, ny);
[xx, yy] = meshgrid(x, y);
rr = sqrt(xx .^ 2 + yy .^ 2);

P = cell(3, 1); P(:) = {zeros(nx, ny)};   % Pattern (P = 1 inside colony; P = 0 ouside colony)
C = cell(3, 1); C(:) = {zeros(nx, ny)};   % Cell density
N = zeros(nx, ny) + N0;                   % Nutrient

ntips0 = ceil(2 * pi * r0 * Density); % initial branch number
ntips0 = max(ntips0, 2);  % initial branch number cannot be less than 2
for j = 1 : 3
    P{j}(rr <= r0) = 1;
    C{j}(P{j} == 1) = C0 * initialFract(j) / (sum(P{j}(:)) * dx * dy);
end
C_pre = C;
ntips = [ntips0, ntips0, ntips0];
Tipx = cell(3, 1); % x coordinates of every tip
Tipy = cell(3, 1); % y coordinates of every tip
Tipx(:) = {zeros(ntips0, totalt / dt_updatebranch + 1)};
Tipy(:) = {zeros(ntips0, totalt / dt_updatebranch + 1)};

dE = zeros(ntips0, 3);
BranchDomain = cell(ntips0, 3); % the domain covered by each branch
BranchDomain(:) = {rr <= r0};

theta = linspace(0, 2 * pi, ntips0 + 1)' + rand * 0;
theta = theta(1 : ntips0);  % growth directions of every branch
theta = theta * ones(1, 3); % one species each column
delta = linspace(-1, 1, 201) * pi;

[MatV1N,MatV2N,MatU1N,MatU2N] = Diffusion(dx,dy,nx,ny,dt,DN); % for solving diffusion equation using ADI method
ib = 0;

for i = 0 : nt

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
        
        ib = ib + 1;
        Tipx{j}(:, ib) = Tipx{j}(:, max(1, ib - 1));
        Tipy{j}(:, ib) = Tipy{j}(:, max(1, ib - 1));
        dBiomass = (C{j} - C_pre{j}) * dx * dy; % total cell growth

        % compute the amount of biomass accumulation in each branch
        BranchDomainSum = cat(3, BranchDomain{:,j});
        BranchDomainSum = sum(BranchDomainSum, 3);
        nn = ntips(j);
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

        jx = j - 1; if jx == 0; jx = 3; end
        [Tipx{j}, Tipy{j}] = tiptracking(Tipx{jx}, Tipy{jx}, Tipx{j}, Tipy{j}, ...
            ib, dl, theta, delta, nn, xx, yy, N, j, noiseamp);
        
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
        ntips(j) = nn;
        if nn > size(Tipx{j}, 1)
            nz = nn - size(Tipx{j}, 1);
            theta = [theta; zeros(nz, 3)];
            BranchDomain = [BranchDomain; cell(nz, 3)];
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