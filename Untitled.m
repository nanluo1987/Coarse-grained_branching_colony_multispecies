
for jj = 1:3

k = 1;
d = sqrt((Tipx{j}(k,ib) - xx) .^ 2 + (Tipy{j}(k,ib) - yy) .^ 2);
ind2 = d <= Width/2;
                
ind = 1 : 2 : nx;
subplot(3, 4, (jj-1)*4+1)
hold off; pcolor(xx(ind, ind), yy(ind, ind), hs{jj}(ind, ind));
shading interp; axis equal;
axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
colorbar
set(gca,'YTick',[], 'XTick',[])
plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
ind = 1 : 2 : nx;

subplot(3, 4, (jj-1)*4+2)
hold off; pcolor(xx(ind, ind), yy(ind, ind), dBiomass{jj}(ind, ind));
shading interp; axis equal;
axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
colorbar
set(gca,'YTick',[], 'XTick',[])
plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
caxis([0 1.8e-4])

subplot(3, 4, (jj-1)*4+3)
hold off; pcolor(xx(ind, ind), yy(ind, ind), dBiomass{jj}(ind, ind).*ind2(ind, ind));
shading interp; axis equal;
axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
colorbar
set(gca,'YTick',[], 'XTick',[])
plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
caxis([0 1.8e-4])

subplot(3, 4, (jj-1)*4+4)
hold off; pcolor(xx(ind, ind), yy(ind, ind), hs{jj}(ind, ind).*dBiomass{jj}(ind, ind).*ind2(ind, ind));
shading interp; axis equal;
axis([-L/2 L/2 -L/2 L/2]); colormap('parula'); hold on
colorbar
set(gca,'YTick',[], 'XTick',[])
plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)

end

figure(1)