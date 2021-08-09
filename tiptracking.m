function [xm2, ym2] = tiptracking(xcell, ycell, ib, dL, theta, delta, nn, xx, yy, N, j, noiseamp)

% input:  coordinates of branch tips
% output: coordinates of branch tips of species j (1 to ib)
% xm1,xm2,ym1,ym2 are matrices (number of tips x number of time points)

xm2 = xcell{j};  
ym2 = ycell{j};

jlist = [1, 2, 3, 1, 2];
xma = xcell{jlist(j + 1)};  xmb = xcell{jlist(j + 2)};
yma = ycell{jlist(j + 1)};  ymb = ycell{jlist(j + 2)};

da = sqrt(diff(xma(:, 1 : ib),1,2) .^ 2 + diff(yma(:, 1 : ib),1,2) .^ 2);
db = sqrt(diff(xmb(:, 1 : ib),1,2) .^ 2 + diff(ymb(:, 1 : ib),1,2) .^ 2);

idx = sum(da, 2) < sum(db, 2);
xm1 = xma; xm1(idx, :) = xmb(idx, :);
ym1 = yma; ym1(idx, :) = ymb(idx, :);

% coordinates of branch tips of species 2 at ib
xv = zeros(size(xm1,1), 1);
yv = zeros(size(xm1,1), 1);

% length of each segments of each branch
d1 = sqrt(diff(xm1(:, 1 : ib),1,2) .^ 2 + diff(ym1(:, 1 : ib),1,2) .^ 2); 
izero = abs(xm1(:, ib)) + abs(ym1(:, ib)) == 0; d1(izero, end) = 0;
d2 = sqrt(diff(xm2(:, 1 : ib - 1),1,2) .^ 2 + diff(ym2(:, 1 : ib - 1),1,2) .^ 2);
L1 = sum(d1, 2);
L2 = sum(d2, 2);

% ind_leader = L2 >= 0;
ilead = L2 + dL  > L1;
ifllw = L2 + dL <= L1;

if ib == 2
    
    xv(ilead) = xm2(ilead, 1) + dL .* sin(theta(ilead,j));
    yv(ilead) = ym2(ilead, 1) + dL .* cos(theta(ilead,j));
    
else
    
    dL0 = dL;
  % ------------------------
  % follow nutrient gradient
    x0 = xm2(:, ib - 1);
    y0 = ym2(:, ib - 1);

    idx = L2 < L1;
    x0(idx) = xm1(idx, ib); x0(idx & izero) = xm1(idx & izero, ib - 1);
    y0(idx) = ym1(idx, ib); y0(idx & izero) = ym1(idx & izero, ib - 1);
    dL(idx) = dL(idx) - (L1(idx) - L2(idx));

    thetaO = ones(nn, 1) * delta;
    TipxO = x0 + dL .* sin(thetaO);
    TipyO = y0 + dL .* cos(thetaO);
    NO = interp2(xx, yy, N, TipxO, TipyO);
    [~, idx] = max(NO, [], 2); % find the direction with maximum nutrient
    TipxO = x0 + dL .* sin(thetaO);
    TipyO = y0 + dL .* cos(thetaO);
    for k = 1 : nn
        if ilead(k)
            xv(k) = TipxO(k, idx(k));
            yv(k) = TipyO(k, idx(k));
            theta(k,j) = thetaO(k, idx(k)) + noiseamp * rand;
        end
    end
    dL = dL0;
    
  % ----------------------------------
  % follow the trajectory of species 1
    cumL1 = cumsum(d1, 2);
    idx = sum(cumL1 < (L2 + dL) * ones(1, size(cumL1, 2)), 2);
    idxlinear = sub2ind(size(d1), 1 : size(d1, 1), idx')';
    dL = dL - (cumL1(idxlinear) - L2);
    
    idx = min(size(d1, 2), idx + 1);
    idxlinear = sub2ind(size(d1), 1 : size(d1, 1), idx')';
    x0 = xm1(idxlinear);
    y0 = ym1(idxlinear);
    
    idx = min(size(d1, 2), idx + 1);
    idxlinear = sub2ind(size(d1), 1 : size(d1, 1), idx')';
    x1 = xm1(idxlinear);
    y1 = ym1(idxlinear);
    
    dx = sqrt((x1 - x0) .^ 2 + (y1 - y0) .^ 2);
    xv(ifllw) = (x1(ifllw) - x0(ifllw)) .* dL(ifllw) ./ dx(ifllw) + x0(ifllw);
    yv(ifllw) = (y1(ifllw) - y0(ifllw)) .* dL(ifllw) ./ dx(ifllw) + y0(ifllw);

end

xm2(:, ib) = xv;
ym2(:, ib) = yv;
    
