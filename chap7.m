%% 7.7
% F = @(t,x) 0;
% ode23tx(F,[0 10],1)
F = @(t,y) [y(2); -y(1)]
ode23tx(F,[0 2*pi],[1;0])