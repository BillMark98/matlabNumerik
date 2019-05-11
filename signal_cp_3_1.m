
% MATLAB script for Computer Problem 3.1.
% Matlab demonstration script for DSB-AM modulation. The message signal % is m(t)=sinc(100t).
echo off
t0=.4;
ts=0.0001;
fc=250;
snr =20;
fs=1/ts;
df=0.3;
t=[-t0:ts:t0];
snr_lin=10^(snr/10);
m=sinc(100*t);
c=cos(2*pi*fc.*t);
u=m.*c;
[M,m,df1]=fftseq(m,ts,df);


M=M/fs;
[U,u,df1]=fftseq(u,ts,df);
U=U/fs;
f=[0:df1:df1*(length(m)- 1)]- fs/2;
% plot the message signal
figure;
plot(t,m(1:length(t)))
xlabel('Time')
% plot the modulated signal.
figure;
plot(t,u(1:length(t)))
xlabel('Time')
% plot the magnitude of the message and the % modulated signal in the frequency domain. figure;
plot(f,abs(fftshift(M)))
xlabel('Frequency')
axis([-1000 1000 0 9*10^(-3)]);
figure;
plot(f,abs(fftshift(U)))
xlabel('Frequency')
axis([-1000 1000 0 4.5*10^(-3)]);