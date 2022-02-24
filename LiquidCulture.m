
C0  = 0.001;   % initial cell density
nt = totalt / dt;  % number of time steps

C = zeros(nt + 1, 3);   % Cell density
N = zeros(nt + 1, 1);   % Nutrient

C(1, :) = C0 * initialFract;
N(1) = N0;  

for i = 1 : nt
    
    fN = N(i) ./ (N(i) + 1) .* (1 - sum(C(i, :)));
    
    aCs = aCs_act;
    if N(i) > N_upper || N(i) < N_lower; aCs = max(aCs_act); end
    
    C(i + 1, :) = C(i, :) + aCs .* fN .* C(i, :) * dt;
    N(i + 1)    = N(i) - bN * fN .* sum(aCs .* C(i, :)) * dt;    

end

Biomass = C(end, :);
BiomassV = C;