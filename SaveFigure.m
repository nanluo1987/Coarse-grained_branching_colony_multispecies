
Ctotal = 5*C{1} + C{2} + C{3};
p1 = 5*C{1}./Ctotal; p1(isnan(p1)) = 0;
p2 = C{2}./Ctotal; p2(isnan(p2)) = 0;
p3 = C{3}./Ctotal; p3(isnan(p3)) = 0;
ind = 1 : 2 : nx;
color1 = [199,61,120]; color2 = [255,192,0]; color3 = [52,117,166];
ColorMap = MarkMixing_3color(color1, color2, color3, p1, p2, p3);

figure(2); clf; set(gcf,'position',[360   227   391   391])
hold off; surf(xx(ind, ind), yy(ind, ind), ones(size(xx(ind, ind))), ColorMap(ind, ind, :))
view([0, 0, 1]); shading interp; axis equal; box on
axis([-L/2 L/2 -L/2 L/2]);
set(gca,'YTick',[], 'XTick',[], 'Visible', 'off')

saveas(gca, [filename '.jpg'])
figure(1)
