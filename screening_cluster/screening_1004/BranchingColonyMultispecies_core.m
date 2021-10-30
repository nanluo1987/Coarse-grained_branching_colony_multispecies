for iter = 1 : 1

speciesName = {'WT','Cheater','Hyperswarmer'}; % name of each species
% other vectors will follow the same order

initialRatio = pyramid_initRatios(iter, :);   % initial ratio of all species
initialFract = initialRatio / sum(initialRatio); % initial fraction of each species

% branch density & width of single-species colonies
Densities = [0.14, 0.14, 0.25];
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
ntips = ntips0;
Tipx = cell(3, 1); % x coordinates of every tip
Tipy = cell(3, 1); % y coordinates of every tip
Tipx(:) = {zeros(ntips0, totalt / dt_updatebranch + 2)};
Tipy(:) = {zeros(ntips0, totalt / dt_updatebranch + 2)};

dE = zeros(ntips0, 3);
BranchDomain = cell(ntips0, 3); % the domain covered by each branch
BranchDomain(:) = {rr <= r0};

theta = linspace(0, 2 * pi, ntips0 + 1)' + rand * 0;
theta = theta(1 : ntips0);  % growth directions of every branch
theta = theta * ones(1, 3); % one species each column
delta = linspace(-1, 1, 201) * pi;

[MatV1N,MatV2N,MatU1N,MatU2N] = Diffusion(dx,dy,nx,ny,dt,DN); % for solving diffusion equation using ADI method
aCs = cell(3, 1);
hs  = cell(3, 1);

for i = 0 : nt
    
    fN = N ./ (N + 1) .* (1 - (C{1} + C{2} + C{3}));
    
    for j = 1 : 3  % j: index of species
        
        % -------------------------------------
        % Cell growth
        aCs{j} = aCs_act(j) * ones(nx, ny); % growth rate as a variable of time and location
        aCs{j}(N > N_upper | N < N_lower) = max(aCs_act); % if nutrient is not within the range grow at full speed
        C{j} = C{j} + aCs{j} .* fN .* C{j} * dt;
    
    end
    
    % -------------------------------------
    % Nutrient distribution
    dN = - bN * fN .* (aCs{1} .* C{1} + aCs{2} .* C{2} + aCs{3} .* C{3}); % Nutrient consumption
    N  = N + dN * dt;
    NV = MatV1N \ (N * MatU1N); N = (MatV2N * NV) / MatU2N; % Nutrient diffusion

    % -------------------------------------
    % Cell allocation and branch extension

    if mod(i, dt_updatebranch / dt) == 0
    
    dBiomass = cell(3, 1);
    for j = find(initialFract > 0) % get growth of each species
        ib = i / (dt_updatebranch / dt) + 2;
        Tipx{j}(:, ib) = Tipx{j}(:, max(1, ib - 1));
        Tipy{j}(:, ib) = Tipy{j}(:, max(1, ib - 1));     
        dBiomass{j} = (C{j} - C_pre{j}) * dx * dy; % total cell growth  
        % gamma (expansion efficiency) = swimming efficiency + swarming efficiency
        % calculate the gamma at each branch tip from the local composition
        nn = ntips;
        tipC = [interp2(xx, yy, C{1}, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib)) ...
                interp2(xx, yy, C{2}, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib)) ...
                interp2(xx, yy, C{3}, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib))];
        tipN =  interp2(xx, yy, N, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib));
        hs{j} = hs_act(j) * ones(nx, ny);
        hs{j}(N > N_upper | N < N_lower) = 0; % swarm only when nutrient is within the range
    end

    for j = find(initialFract > 0) % get pattern of each species
        
        % update width and density with time
        currFracts = [sum(C{1}, 'all') sum(C{2}, 'all') sum(C{3}, 'all')];
        currFracts = currFracts ./ sum(currFracts); %single vector, the below is spatial grid. 
        Density = currFracts * Densities';
        Width   = currFracts * Widths';      
        
        for jj = find(initialFract > 0) % get growth of each species at the tip 
            for k = 1 : nn
                d = sqrt((Tipx{j}(k,ib) - xx) .^ 2 + (Tipy{j}(k,ib) - yy) .^ 2);
                ind = d <= Width/2;
                dE(k,jj) = sum(hs{jj} .* dBiomass{jj} .* ind,'all'); % the growth of the jjth species within the kth branch of the jth species
            end
        end
        
        % extension rate of each branch  
        aCs_tip = interp2(xx, yy, aCs{j}, Tipx{j}(1:nn,ib), Tipy{j}(1:nn,ib));
        dl = (gs(j) * 1 + aCs_tip .* sum(dE(1:nn,:), 2))./ Width;
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
                TipxNew(nn,ib) = Tipx{j}(k,ib) + 0.7654*dl(k) * sin(theta(k,j) + 0.625 * pi); % splitting the old tip to two new tips
                TipyNew(nn,ib) = Tipy{j}(k,ib) + 0.7654*dl(k) * cos(theta(k,j) + 0.625 * pi); % numbers are to ensure extension = dl & separation angle = 90 after splitting
                TipxNew(k,ib) = TipxNew(k,ib) + 0.7654*dl(k) * sin(theta(k,j) - 0.625 * pi);
                TipyNew(k,ib) = TipyNew(k,ib) + 0.7654*dl(k) * cos(theta(k,j) - 0.625 * pi);
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

        idx = abs(Tipx{j}(:, ib)) > L/2 | abs(Tipy{j}(:, ib)) > L/2;
        Tipx{j}(idx, ib) = Tipx{j}(idx, ib - 1);
        Tipy{j}(idx, ib) = Tipy{j}(idx, ib - 1);

    %         %     Growth stops when approaching edges
    %         TipR = sqrt(Tipx{j}(:,ib).^2 + Tipy{j}(:,ib).^2);
    %         idx = TipR > 0.9*L/2;
    %         Tipx{j}(idx, ib) = Tipx{j}(idx, ib - 1);
    %         Tipy{j}(idx, ib) = Tipy{j}(idx, ib - 1);

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
        C_relo = frac_relo * sum(dBiomass{j}(:)) / sum(Capacity(:)) * Capacity / (dx * dy);
        C{j} = C_pre{j} + C_relo + (1 - frac_relo) * dBiomass{j} / (dx * dy);
        
    end
    
    C_pre = C;

    end
    

end

Output_Biomass((iter - 1) * 3 + 1 : iter * 3) = [sum(C{1},'all'); sum(C{2},'all'); sum(C{3},'all')];
Sizes = zeros(3, 1);
for j = 1 : 3
    if ~isempty(rr(C{j} > 0))
        Sizes(j) = max(rr(C{j} > 0),[],'all');
    else
        Sizes(j) = 0;
    end
end
Output_Sizes((iter - 1) * 3 + 1 : iter * 3) = Sizes;
fprintf('iter = %d\n', iter)

end