%% 6.9
f = inline('cos(x).^2')
quadtx(f,0,4*pi)
2 * pi
quad(f,0,4 * pi)
ff = @(x) cos(x).^2
integral(ff,0,4 * pi)