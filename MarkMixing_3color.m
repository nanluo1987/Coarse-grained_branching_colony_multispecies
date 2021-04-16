function colors = MarkMixing_3color(color1, color2, color3, cx, cy, cz)
% Generating a smooth color map using three colors
% Based on Mark's method for generating a color gradient between two colors
% https://stackoverflow.com/questions/22607043/color-gradient-algorithm

gamma = .43;
color1_lin = from_sRGB(color1);
bright1 = sum(color1_lin)^gamma;
color2_lin = from_sRGB(color2);
bright2 = sum(color2_lin)^gamma;
color3_lin = from_sRGB(color3);
bright3 = sum(color3_lin)^gamma;
colors = zeros(size(cx,1), size(cx,2), 3);
intensity = lerp(bright1, bright2, bright3, cx, cy, cz) .^ (1/gamma);
for dim = 1 : 3    
    color = lerp(color1_lin(dim), color2_lin(dim), color3_lin(dim), cx, cy, cz);
    colors(:, :, dim) = color;
end
sumcolor = colors(:, :, 1) + colors(:, :, 2) + colors(:, :, 3);
sumcolor(sumcolor == 0) = intensity(sumcolor == 0);
for dim = 1 : 3
    colors(:, :, dim) = colors(:, :, dim) .* intensity ./ sumcolor;
    colors(:, :, dim) = to_sRGB_f(colors(:, :, dim));
end

function f = to_sRGB_f(x)
%     ''' Returns a sRGB value in the range [0,1]
%         for linear input in [0,1].
f = 12.92*x;
f(x > 0.0031308) = (1.055 * (x(x > 0.0031308) .^ (1/2.4))) - 0.055;

function f = to_sRGB(x)
%     ''' Returns a sRGB value in the range [0,255]
%         for linear input in [0,1]
f = round(255.9999 * to_sRGB_f(x));

function y = from_sRGB(x)
%     ''' Returns a linear value in the range [0,1]
%         for sRGB input in [0,255].
x = x / 255.0;
y = x / 12.92;
y(x > 0.04045) = ((x(x > 0.04045) + 0.055) / 1.055) .^ 2.4;

function f = lerp(color1, color2, color3, cx, cy, cz)
f = color1 * cx + color2 * cy + color3 * cz;

