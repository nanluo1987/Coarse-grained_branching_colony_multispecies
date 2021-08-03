function [xm2, ym2] = tiptracking(xm1, ym1, xm2, ym2, ib, dl, theta, delta, nn, xx, yy, N, j, noiseamp)

% input:  coordinates of branch tips of species 1 (1 to ib) & 2 (1 to ib - 1)
%         dl of species 2
% output: coordinates of branch tips of species 2 (1 to ib)
% xm1,xm2,ym1,ym2 are matrices (number of tips x number of time points)

% coordinates of branch tips of species 2 at ib
xv = zeros(size(xm1,1), 1);
yv = zeros(size(xm1,1), 1);

% length of each segments of each branch
d1 = sqrt(diff(xm1(:, 1 : ib),1,2) .^ 2 + diff(ym1(:, 1 : ib),1,2) .^ 2);
d2 = sqrt(diff(xm2(:, 1 : ib - 1),1,2) .^ 2 + diff(ym2(:, 1 : ib - 1),1,2) .^ 2);
L1 = sum(d1, 2);
L2 = sum(d2, 2);

ind_leader = L2 >= 0;
% ind_leader = L2 + dl  > L1;
% ind_fllwer = L2 + dl <= L1;

if ib == 1
    xv(ind_leader) = xm2(ind_leader, 1) + dl .* sin(theta(ind_leader,j));
    yv(ind_leader) = ym2(ind_leader, 1) + dl .* cos(theta(ind_leader,j));
else
    
    % ------------------------
    % follow nutrient gradient
    x0 = xm2(:, ib - 1);
    y0 = ym2(:, ib - 1);

    ind = L2 < L1;
    x0(ind) = xm1(ind, ib);
    y0(ind) = ym1(ind, ib);
    dl(ind) = dl(ind) - (L2(ind) - L1(ind));


        thetaO = ones(sum(ind_leader), 1) * delta;
        TipxO = x0(ind_leader) + dl(ind_leader) .* sin(thetaO);
        TipyO = y0(ind_leader) + dl(ind_leader) .* cos(thetaO);
        NO = interp2(xx, yy, N, TipxO, TipyO);
        [~, ind] = max(NO, [], 2); % find the direction with maximum nutrient
        TipxO = x0(ind_leader) + dl(ind_leader) .* sin(thetaO);
        TipyO = y0(ind_leader) + dl(ind_leader) .* cos(thetaO);
        for k = 1 : nn
            if ind_leader(k)
                xv(k) = TipxO(k, ind(k));
                yv(k) = TipyO(k, ind(k));
                theta(k,j) = thetaO(k, ind(k)) + noiseamp * rand;
            end
        end

    % ----------------------------------
    % follow the trajectory of species 1

end


xm2(:, ib) = xv;
ym2(:, ib) = yv;
    
