omega = 1;
beta = 1.3;
gamma = 0.2;
z = linspace(-10,0,100);
figure(1)
for t = 0:0.1 :20
y1 = cos(omega*t - beta * z);
y2 = gamma * cos(omega*t + beta * z);
y3 = y1 + y2;
subplot(3,1,1);
plot(z,y1);
ylim([-2,2]);
xlim([-10,0]);

subplot(3,1,2);
plot(z,y2);
ylim([-2,2]);
xlim([-10,0]);

subplot(3,1,3);
plot(z,y3);
ylim([-2,2]);
xlim([-10,0]);

pause(0.001);
end