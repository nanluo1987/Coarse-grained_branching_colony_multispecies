clear
figure(1)
taskID = str2num(getenv('SLURM_ARRAY_TASK_ID'));
for iter = 1:7
%% Parameters
L      = 90;    % domain size
totalt = 12;    % total time
dt     = 0.02;  % time step
nx     = 1001; ny = nx; % number of

speciesName = {'WT','Cheater','Hyperswarmer'}; % name of each species
% other vectors will follow the same order

pyramid_initRatios = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1];
initialRatio = pyramid_initRatios(iter, :);   % initial ratio of all species
initialFract = initialRatio / sum(initialRatio); % initial fraction of each species
filename = "gamx_" + string(taskID) + "_" + strjoin(string(initialRatio), "") + '.jpg';


aCs = [1, 1.1, 1] * 1.5;      % cell growth rate of each species
gs = 2*[1 1 2]*2;                 % swimming motility
h1s = 22 *2;             % swarming motility coefficient of WT
h3s = 20 *2;             % swarming motility coefficient of hyperswarmer

bN = 150;   % nutrient consumption rate
DN = 7;     % nutrient diffusivity
KN = 1.2;   % half-saturation conc of nutrient-dependent growth
N0 = 8;     % initial nutrient conc.
Cmax = 0.2; % cell carrying capacity
noiseamp = 0; % noise amplitude of branch direction

% branch density & width of single-species colonies
Densities = [0.14, 0.14, 0.2];
Widths    = [5, 5, 20];

% branch density & width of mixed colonies
Density = initialFract * Densities';
Width   = initialFract * Widths';

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

r0  = 5;     % initial radius
C0  = 1.6;   % initial cell density
ntips0 = ceil(2 * pi * r0 * Density); % initial branch number
ntips0 = max(ntips0, 2);  % initial branch number cannot be less than 2
for j = 1 : 3
    P{j}(rr <= r0) = 1;
    C{j}(P{j} == 1) = C0 * initialFract(j) / (sum(P{j}(:)) * dx * dy);
end
C_pre = C;
ntips = [ntips0, ntips0, ntips0];
Tipx = zeros(ntips0, 3); % x coordinates of every tip
Tipy = zeros(ntips0, 3); % y coordinates of every tip

dE = zeros(ntips0, 3);
BranchDomain = cell(ntips0, 3); % the domain covered by each branch
BranchDomain(:) = {rr <= r0};

theta = linspace(0, 2 * pi, ntips0 + 1)' + rand * 0;
theta = theta(1 : ntips0);  % growth directions of every branch
theta = theta * ones(1, 3); % one species each column
delta = linspace(-1, 1, 201) * pi;

[MatV1N,MatV2N,MatU1N,MatU2N] = Diffusion(dx,dy,nx,ny,dt,DN); % for solving diffusion equation using ADI method

for i = 0 : nt

    % -------------------------------------
    % Nutrient distribution

    fN = N ./ (N + KN) .* (1 - (C{1} + C{2} + C{3}) / Cmax);
    dN = - bN * fN .* (aCs(1) * C{1} + aCs(2) * C{2} + aCs(3) * C{3}); % Nutrient consumption
    N  = N + dN * dt;
    NV = MatV1N \ (N * MatU1N); N = (MatV2N * NV) / MatU2N; % Nutrient diffusion

    for j = [1, 3, 2]  % j: index of species (calculate cheater lastly to ensure cheater is not the fastest one)

    % -------------------------------------
    % Cell growth

    C{j} = C{j} + aCs(j) * fN .* C{j} * dt;

    % -------------------------------------
    % Cell allocation and branch extension

    if mod(i, 0.1/dt) == 0 && initialFract(j) > 0

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

        % gamma (expansion efficiency) = swimming efficiency + swarming efficiency
%         currFract = cellfun(@(x) sum(x, 'all'), C);
%         currFract = currFract ./ sum(currFract);
%         gammas = [gs(1)+h1s(1)*currFract(1)+h3s(1)*currFract(3), ...
%             gs(2)+h1s(2)*currFract(1)+h3s(2)*currFract(3), ...
%             gs(3)+h1s(3)*currFract(1)+h3s(3)*currFract(3)];
        currFract_wt = C{1} ./ (C{1} + C{2} + C{3});
        currFract_hs = C{3} ./ (C{1} + C{2} + C{3});
        tipFract_wt = interp2(xx, yy, currFract_wt, Tipx(1:nn,j), Tipy(1:nn,j));
        tipFract_hs = interp2(xx, yy, currFract_hs, Tipx(1:nn,j), Tipy(1:nn,j));
        gammas = gs(j) + h1s * tipFract_wt + h3s * tipFract_hs;
        
        % extension rate of each branch  
%         dl = gammas(j) * dE(1:nn,j) ./ Width;
        dl = gammas .* dE(1:nn,j) ./ Width;
        if i == 0; dl = 0.5; end

        % Bifurcation
        R = 3/2 / Density;  % a branch will bifurcate if there is no other branch tips within the radius of R
        TipxNew = Tipx(:,j); TipyNew = Tipy(:,j);
        thetaNew = theta(:,j); dlNew = dl;
        BranchDomainNew = BranchDomain(:,j);
        for k = 1 : nn
            dist2othertips = sqrt((TipxNew - Tipx(k,j)) .^ 2 + (TipyNew - Tipy(k,j)) .^ 2);
            dist2othertips = sort(dist2othertips);
            if dist2othertips(2) > R
                nn = nn + 1;
                TipxNew(nn) = Tipx(k,j) + dl(k) * sin(theta(k,j) + 0.5 * pi); % splitting the old tip to two new tips
                TipyNew(nn) = Tipy(k,j) + dl(k) * cos(theta(k,j) + 0.5 * pi);
                TipxNew(k) = TipxNew(k) + dl(k) * sin(theta(k,j) - 0.5 * pi);
                TipyNew(k) = TipyNew(k) + dl(k) * cos(theta(k,j) - 0.5 * pi);
                dlNew(nn) = dl(k) / 2;
                dlNew(k)  = dl(k) / 2;
                thetaNew(nn) = theta(k,j);
                BranchDomainNew{nn} = BranchDomain{k,j};
            end
        end
        ntips(j) = nn;
        if nn > size(Tipx, 1)
            nz = nn - size(Tipx, 1);
            Tipx = [Tipx; zeros(nz, 3)];
            Tipy = [Tipy; zeros(nz, 3)];
            theta = [theta; zeros(nz, 3)];
            BranchDomain = [BranchDomain; cell(nz, 3)];
        end
        Tipx(:,j) = TipxNew; Tipy(:,j) = TipyNew;
        theta(:,j) = thetaNew; dl = dlNew;
        BranchDomain(:,j) = BranchDomainNew;      
        
        % Determine branch extension directions
        if i == 0
            Tipx(1:nn,j) = Tipx(1:nn,j) + dl .* sin(theta(1:nn,j));
            Tipy(1:nn,j) = Tipy(1:nn,j) + dl .* cos(theta(1:nn,j));
        else
            thetaO = ones(nn, 1) * delta;
            TipxO = Tipx(1:nn,j) + dl .* sin(thetaO);
            TipyO = Tipy(1:nn,j) + dl .* cos(thetaO);
            NO = interp2(xx, yy, N, TipxO, TipyO);
            [~, ind] = max(NO, [], 2); % find the direction with maximum nutrient
            TipxO = Tipx(1:nn,j) + dl .* sin(thetaO);
            TipyO = Tipy(1:nn,j) + dl .* cos(thetaO);
            for k = 1 : nn
                Tipx(k,j) = TipxO(k, ind(k));
                Tipy(k,j) = TipyO(k, ind(k));
                theta(k,j) = thetaO(k, ind(k)) + noiseamp * rand;
            end
        end
        
        % Growth stops when approaching edges
%         ind = TipR > 0.85 * L/2;
%         Tipx(ind) = Tipx_pre(ind);
%         Tipy(ind) = Tipy_pre(ind);

        % Fill the width of the branches
        for k = 1 : nn
            d = sqrt((Tipx(k,j) - xx) .^ 2 + (Tipy(k,j) - yy) .^ 2);
            P{j}(d <= Width/2) = 1;
            BranchDomain{k,j} = BranchDomain{k,j} | (d <= Width/2);
        end

        % relocate dBiomass
        Capacity = Cmax - C_pre{1} - C_pre{2} - C_pre{3};    % remaining cell capacity
        Capacity(P{j} == 0) = 0;   % no capacity outside the colony
        C_relo = sum(dBiomass(:)) / sum(Capacity(:)) * Capacity / (dx * dy);
        C{j} = C_pre{j} + C_relo;
        C_pre{j} = C{j};
    end 
    end
end  
% Plot each species
ind = 1 : 2 : nx;
for j = 1:3
subplot(2, 3, j)
    hold off; pcolor(xx(ind, ind), yy(ind, ind), C{j}(ind, ind));
    shading interp; axis equal;
    axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
    colorbar
    set(gca,'YTick',[], 'XTick',[])
    plot(Tipx(:,j), Tipy(:,j), '.', 'markersize', 5)
    title(speciesName{j})
    drawnow
end
%Plot all species & save. 
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
    plot(Tipx(:,j), Tipy(:,j), '.', 'markersize', 5)
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
saveas(gca,  "1_" + filename);

% save results 
figure(2); clf
Ctotal = C{1} + C{2} + C{3};
p1 = C{1}./Ctotal; p1(isnan(p1)) = 0;
p2 = C{2}./Ctotal; p2(isnan(p2)) = 0;
p3 = C{3}./Ctotal; p3(isnan(p3)) = 0;
ind = 1 : 2 : nx;
color1 = [199,61,120]; color2 = [255,192,0]; color3 = [52,117,166];
ColorMap = MarkMixing_3color(color1, color2, color3, p1, p2, p3);
hold off; surf(xx(ind, ind), yy(ind, ind), ones(size(xx(ind, ind))), ColorMap(ind, ind, :))
view([0, 0, 1]); shading interp; axis equal; box on
axis([-L/2 L/2 -L/2 L/2]);
set(gca,'YTick',[], 'XTick',[])
saveas(gca,  "2_" + filename);

end