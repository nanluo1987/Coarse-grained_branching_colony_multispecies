clear

%% Parameters
L      = 90;    % domain size
totalt = 16;    % total time
dt     = 0.02;  % time step
nx     = 1001; ny = nx; % number of nodes

speciesName = {'WT','Cheater','Hyperswarmer'}; % name of each species
% other vectors will follow the same order

initialRatio = [1, 1, 1];   % initial ratio of all species
initialFract = initialRatio / sum(initialRatio); % initial fraction of each species

bN = 25;             % nutrient consumption rate
DN = 7;              % nutrient diffusivity
N0 = 8/1.2;         % initial nutrient conc.

aCs = [1, 1.1, 0.9] * 1.5/2;      % cell growth rate of each species
gs = 1*[1 1 2];               % swimming motility
hs = [1, 0, 0.9] * 4;        % swarming motility coefficients

% branch density & width of single-species colonies
Densities = [0.14, 0.14, 0.2];
Widths    = [5, 5, 20];

% branch density & width of mixed colonies
Density = initialFract * Densities';
Width   = initialFract * Widths';

r0  = 5;   % initial radius
C0  = 8;   % initial cell density

noiseamp = 0;        % noise amplitude of branch direction

%% Initialization

nt = totalt / dt;  % number of time steps
dx = L / (nx - 1); dy = dx;
x  = linspace(-L/2, L/2, nx);
y  = linspace(-L/2, L/2, ny);
[xx, yy] = meshgrid(x, y);
rr = sqrt(xx .^ 2 + yy .^ 2);

P = zeros(nx, ny);   % Pattern (P = 1 inside colony; P = 0 ouside colony)
C = cell(3, 1); C(:) = {zeros(nx, ny)};   % Cell density
N = zeros(nx, ny) + N0;                   % Nutrient
Rad = zeros(3, 1) + r0;     % Radius of each species

ntips0 = ceil(2 * pi * r0 * Density); % initial branch number
ntips0 = max(ntips0, 2);  % initial branch number cannot be less than 2
for j = 1 : 3
    P(rr <= r0) = 1;
    C{j}(P == 1) = C0 * initialFract(j) / (sum(P(:)) * dx * dy);
end
C_pre = C;
ntips = ntips0;
Tipx = zeros(ntips0, 1); % x coordinates of every tip
Tipy = zeros(ntips0, 1); % y coordinates of every tip

theta = linspace(0, 2 * pi, ntips0 + 1)' + rand * 0;
theta = theta(1 : ntips0);  % growth directions of every branch
delta = linspace(-1, 1, 201) * pi;

Tipx = r0 * sin(theta);
Tipy = r0 * cos(theta);

[MatV1N,MatV2N,MatU1N,MatU2N] = Diffusion(dx,dy,nx,ny,dt,DN); % for solving diffusion equation using ADI method
C_front = zeros(3,1);
Rad_pre = Rad;

for i = 0 : nt

    % -------------------------------------
    % Nutrient distribution

    fN = N ./ (N + 1) .* (1 - (C{1} + C{2} + C{3}));
    dN = - bN * fN .* (aCs(1) * C{1} + aCs(2) * C{2} + aCs(3) * C{3}); % Nutrient consumption
    N  = N + dN * dt;
    NV = MatV1N \ (N * MatU1N); N = (MatV2N * NV) / MatU2N; % Nutrient diffusion
    
    % gamma (expansion efficiency) = swimming efficiency + swarming efficiency
    % calculate the local gamma at the front of each species
    for j = 1 : 3
        C_front(j) = sum(sum(C{j}(rr <= Rad(j) & rr > Rad(j) - 2 & P == 1)));
    end
    frontFract = C_front ./ sum(C_front);
    gammas = gs + sum(hs .* frontFract', 2);
    
    for j = 1 : 3  % j: index of species

        % -------------------------------------
        % Cell growth
        dC = aCs(j) * fN .* C{j} * dt;
        C{j} = C{j} + dC;
    
    end
    
    % -------------------------------------
    % Branch extension and biomass relocation

    if mod(i, 0.2/dt) == 0
        
        % Update width and density with time
        currFracts = [sum(C{1}, 'all') sum(C{2}, 'all') sum(C{3}, 'all')];
        currFracts = currFracts ./ sum(currFracts); %single vector, the below is spatial grid. 
        Density = currFracts * Densities';
        Width   = currFracts * Widths';
        
        for j = 1 : 3
            % Front advance
            frontRingx = Rad(j) * sin(delta);
            frontRingy = Rad(j) * cos(delta);
            frontP = interp2(xx, yy, P, frontRingx, frontRingy);
            totalWidth = sum(frontP > 0.5) / sum(frontP) * 2 * pi * Rad(j);
            dr = gammas(j) * sum(C{j} - C_pre{j},'all') * dx * dy / totalWidth;
            %         if i == 0; dr = 0.5; end
            Rad(j) = Rad(j) + dr;
        end
        
        % branch extension dl = max(dr) front advance rate of the leading species
        dl = max(Rad - Rad_pre);
        Rad_pre = Rad;
        
        % Bifurcation
        R = 3/2 / Density;  % a branch will bifurcate if there is no other branch tips within the radius of R
        TipxNew = Tipx; TipyNew = Tipy;
        thetaNew = theta;
        nn = length(Tipx);
        for k = 1 : nn
            dist2othertips = sqrt((TipxNew - Tipx(k)) .^ 2 + (TipyNew - Tipy(k)) .^ 2);
            dist2othertips = sort(dist2othertips);
            if dist2othertips(2) > R
                nn = nn + 1;
                TipxNew(nn) = Tipx(k) + dl * sin(theta(k) + 0.5 * pi); % splitting the old tip to two new tips
                TipyNew(nn) = Tipy(k) + dl * cos(theta(k) + 0.5 * pi);
                TipxNew(k) = TipxNew(k) + dl * sin(theta(k) - 0.5 * pi);
                TipyNew(k) = TipyNew(k) + dl * cos(theta(k) - 0.5 * pi);
                thetaNew(nn) = theta(k);
            end
        end
        Tipx = TipxNew; Tipy = TipyNew; theta = thetaNew;
        
        % Determine branch extension directions
        if i == 0
            Tipx = Tipx + dl .* sin(theta);
            Tipy = Tipy + dl .* cos(theta);
        else
            thetaO = ones(nn, 1) * delta;
            TipxO = Tipx + dl .* sin(thetaO);
            TipyO = Tipy + dl .* cos(thetaO);
            NO = interp2(xx, yy, N, TipxO, TipyO);
            [~, ind] = max(NO, [], 2); % find the direction with maximum nutrient
            TipxO = Tipx + dl .* sin(thetaO);
            TipyO = Tipy + dl .* cos(thetaO);
            for k = 1 : nn
                Tipx(k) = TipxO(k, ind(k));
                Tipy(k) = TipyO(k, ind(k));
                theta(k) = thetaO(k, ind(k)) + noiseamp * rand;
            end
        end

        % Fill the width of the branches
        for k = 1 : nn
            d = sqrt((Tipx(k) - xx) .^ 2 + (Tipy(k) - yy) .^ 2);
            P(d <= Width/2) = 1;
        end

        % Relocate biomass
        for j = 1 : 3
            if Rad(j) < max(Rad)
                Pj = P == 1 & rr <= Rad(j) | rr <= r0;
            else
                Pj = P == 1;
            end
            Capacity = 1 - C_pre{1} - C_pre{2} - C_pre{3};    % remaining cell capacity
            Capacity(Pj == 0) = 0;   % no capacity outside the colony
            frac_relo = 1;
            C_relo = frac_relo * sum(C{j} - C_pre{j}, 'all') / sum(Capacity(:)) * Capacity; % !!!no dx*dy here
            C{j} = C_pre{j} + C_relo + (1 - frac_relo) * sum(C{j} - C_pre{j}, 'all');
            C_pre{j} = C{j};
            
            % Plot each species
            ind = 1 : 2 : nx;
            subplot(2, 3, j)
                hold off; pcolor(xx(ind, ind), yy(ind, ind), C{j}(ind, ind));
                shading interp; axis equal;
                axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
                colorbar
                set(gca,'YTick',[], 'XTick',[])
                plot(Tipx, Tipy, '.', 'markersize', 5)
                title(speciesName{j})
                drawnow
        end
            
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
            plot(Tipx, Tipy, '.', 'markersize', 5)
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
    
    % Growth stops when approaching edges
%     TipR = sqrt(Tipx.^2 + Tipy.^2);
%     if max(TipR(:)) > 0.9 * L/2
%         break
%     end