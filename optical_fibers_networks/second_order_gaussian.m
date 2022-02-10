%% This Matlab code aims to produce simulations of signals transmitting ...
% into single mode optical fibers. The first part is dedicated to
% propagation of gaussian pulses while the second part simulates the
% propagation of random bit sequence.

%% 1.1.a-second order dispersion only, gaussian shape, no loss, analytical formula
close all;

% constants and parameters all units in SI units

D = (20*10^-12)/(10^-9*10^3); %dispersion parameter : s/m^2
S = (0.07*10^-12)/(10^-18*10^3); %dispersion slope : s/m^3
lambda0 = 1500*10^-9; %wavelength : m
c = 299792458; %speed of light : m/s
f0 = c/lambda0;
w0 = 2*pi*f0;
pw = 10*10^-12; %pulse width : s
length = 100*1000; %length of fiber : m
beta2 = -D*lambda0^2/(2*pi*c); %propagation constant of second order
beta3 = S*(lambda0^2/(2*pi*c))^2; %propagation constant of third order
T0 = pw/1.665; % pulse width T0

% The field at a distance z, as a function of the time t is given by the
% integral computed below

wind = 300;
n = 1001;
t = linspace(-wind*pw,wind*pw,n); % time values on which we compute the gaussian
Ez = zeros(1,n);
int1 = zeros(1,n);
it = f0; %boundaries of the integral
k = 2*T0^2;

%initial field
Ei = exp(-t.^2/(2*T0^2));

%we compute the output field for second order dispersion here the Fourier
%transform of the input was computed analytically
fun1 = @(dw,c) sqrt(k*pi)*exp(-pi^2*k*(dw.^2)).*exp((1i/2)*beta2*length*dw.^2).*exp(-1i*dw*c);
for i = 1:n
    int1(i) = integral(@(dw)fun1(dw,t(i)),-inf,inf);
end
Ez1 = (1/(2*pi))*int1;
% Ez1 = abs(Ez1);

figure;
subplot(2,2,1);
plot(t,Ei);
title('Input signal into the fiber');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Ein/Einmax[V/m]');

% figure;
% plot(t,abs(int1));
% title('Integral of the signal out of the fiber for second order');
% xlabel('time [s]');
% ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

% figure;
subplot(2,2,2);
plot(t,abs(Ez1));
title('Output signal out of the fiber for second order - Gaussian - Analytic');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

%% 1.1.b-second order dispersion only, gaussian shape, no loss, numerical

%we also compute the integration using fft, it is a numerical resolution
ns = 24000; %number of sample
ns = 2^nextpow2(ns); %512
ws = 300;

%temporal domain operations
ts = linspace(-ws*pw,ws*pw,ns); % time values on which we compute the gaussian
Eis = exp(-ts.^2/(2*T0^2)); %gaussian values
Ezs1 = zeros(1,ns); %pre allocation of the output field

%frequency domain operation

fs = ns/(2*ws*pw); % fsample : sample at (number of indices/total window duration)
delta = 1/(2*ws*pw); %interval between fft frequencies
fshift = (-ns/2:ns/2-1)*(fs/ns); % zero-centered frequency range

% To use the fft function to convert the signal to the frequency domain,
% first identify a new input length that is the next power of 2 from the 
% original signal length. This will pad the signal X with trailing zeros 
% in order to improve the performance of fft.
% np = 2^nextpow2(ns); 
% 
% fs = np/(2*ws*pw); % fsample : sample at (number of indices/total window duration)
% delta = 1/(2*ws*pw); %interval between fft frequencies
% fshift = (-np/2:np/2-1)*(fs/np); % zero-centered frequency range

Fi = fftshift(fft(Eis,ns));
Fi = abs(Fi/ns);

% figure;
subplot(2,2,3);
plot(fshift,Fi);
title('Fourier Transform of the input signal to the fiber for second order- Gaussian - Numeric');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

%integral using fft

% V1 = exp((1i/2)*beta2*length*fshift.^2).*exp(-1i*fshift*c);
for i = 1:ns
   Ezs1(i) = (1/(2*pi))*2*pi*delta*sum(Fi.*exp((1i/2)*beta2*length*fshift.^2).*exp(-1i*fshift*ts(i)));
end
Ezs1 =  Ezs1/(delta*(2*pi));

% figure;
subplot(2,2,4);
plot(ts,abs(Ezs1));
title('Output signal out of the fiber for second order - Gaussian - Numerical');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

%% 1.2.a-third order dispersion only, gaussian shape, no loss, analytical

wind2 = 300;
t2 = linspace(-wind2*pw,wind2*pw,n);
int2 = zeros(1,n);

%we compute the output field for second order dispersion
fun2 = @(dw,c) sqrt(k*pi)*exp(-pi^2*k*(dw.^2)).*exp((1i/6)*beta3*length*dw.^3).*exp(-1i*dw*c);
for i = 1:n
    int2(i) = integral(@(dw)fun2(dw,t2(i)),-inf,inf);
end
Ez2 = (1/(2*pi))*int2;
%Ez2 = abs(Ez2);

figure;
subplot(2,2,1);
plot(t2,Ei);
title('Input signal into the fiber');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Ein/Einmax[V/m]');

% figure;
% plot(t2,abs(int2));
% title('Integral of the signal out of the fiber for third order');
% xlabel('time [s]');
% ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

% figure;
subplot(2,2,2);
plot(t2,abs(Ez2));
title('Output signal out of the fiber for third order- Gaussian - Analytic');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

% Diff = abs(Ez1-Ez2);
% plot(t,Diff)
%% 1.2.b-third order dispersion only, gaussian shape, no loss, numerical

%we also compute the integration using fft, it is a numerical resolution
ws2 = 300;

%temporal domain operations
ts2 = linspace(-ws2*pw,ws2*pw,ns); % time values on which we compute the gaussian
Eis2 = exp(-ts2.^2/(2*T0^2)); %gaussian values
Ezs2 = zeros(1,ns); %pre allocation of the output field

%frequency domain operation

fs2 = ns/(2*ws2*pw); % fsample : sample at (number of indices/total window duration)
delta2 = 1/(2*ws2*pw); %interval between fft frequencies
fshift2 = (-ns/2:ns/2-1)*(fs2/ns); % zero-centered frequency range

% To use the fft function to convert the signal to the frequency domain,
% first identify a new input length that is the next power of 2 from the 
% original signal length. This will pad the signal X with trailing zeros 
% in order to improve the performance of fft.
% np = 2^nextpow2(ns); 
% 
% fs = np/(2*ws*pw); % fsample : sample at (number of indices/total window duration)
% delta = 1/(2*ws*pw); %interval between fft frequencies
% fshift = (-np/2:np/2-1)*(fs/np); % zero-centered frequency range

Fi2 = fftshift(fft(Eis2,ns));
Fi2 = abs(Fi2/ns);

% figure;
subplot(2,2,3);
plot(fshift2,Fi2);
title('Fourier Transform of the input signal to the fiber for third order - Gaussian - Numeric');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

%integral using fft

% V1 = exp((1i/2)*beta2*length*fshift2.^2).*exp(-1i*fshift2*c);
for i = 1:ns
   Ezs2(i) = (1/(2*pi))*2*pi*delta2*sum(Fi2.*exp((1i/6)*beta3*length*fshift2.^3).*exp(-1i*fshift2*ts2(i)));
end
Ezs2 =  Ezs2/(delta2*(2*pi));

% figure;
subplot(2,2,4);
plot(ts2,abs(Ezs2));
title('Output signal out of the fiber for third order- Gaussian - Numeric');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

% figure;
% Diff1 = abs(Ez1-Ez2);
% plot(t,Diff1)
% 
% figure
% Diff2 = abs(Ezs1-Ezs2);
% plot(ts,Diff2)
%% 2.1.a-second order dispersion only, super-gaussian shape, no loss, numerical

%we also compute the integration using fft, it is a numerical resolution
ws3 = 300;

%temporal domain operations
ts3 = linspace(-ws3*pw,ws3*pw,ns); % time values on which we compute the gaussian
% Eis3 = exp(-ts3.^2/(2*T0^2)); %gaussian values
Eis3 = exp(-0.5*(ts3/T0).^6); %super-gaussian input values
Ezs3 = zeros(1,ns); %pre allocation of the output field

figure;
subplot(3,2,[1,2]);
plot(ts3,Eis3);
title('Input signal to the fiber for Super Gaussian');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

%frequency domain operation

fs3 = ns/(2*ws3*pw); % fsample : sample at (number of indices/total window duration)
delta3 = 1/(2*ws3*pw); %interval between fft frequencies
fshift3 = (-ns/2:ns/2-1)*(fs3/ns); % zero-centered frequency range

% To use the fft function to convert the signal to the frequency domain,
% first identify a new input length that is the next power of 2 from the 
% original signal length. This will pad the signal X with trailing zeros 
% in order to improve the performance of fft.
% np = 2^nextpow2(ns); 
% 
% fs = np/(2*ws*pw); % fsample : sample at (number of indices/total window duration)
% delta = 1/(2*ws*pw); %interval between fft frequencies
% fshift = (-np/2:np/2-1)*(fs/np); % zero-centered frequency range

Fi3 = fftshift(fft(Eis3,ns));
Fi3 = abs(Fi3/ns);

% figure;
subplot(3,2,3);
plot(fshift3,Fi3);
title('Fourier Transform of the Super Gaussian to the fiber for second order');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

%integral using fft

% V1 = exp((1i/2)*beta2*length*fshift2.^2).*exp(-1i*fshift2*c);
for i = 1:ns
   Ezs3(i) = (1/(2*pi))*2*pi*delta3*sum(Fi3.*exp((1i/2)*beta2*length*fshift3.^2).*exp(-1i*fshift3*ts3(i)));
end
Ezs3 =  Ezs3/(delta3*(2*pi));

% figure;
subplot(3,2,4);
plot(ts3,abs(Ezs3));
title('Output signal out of the fiber for second order for Super Gaussian');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

%% 2.1.a-third order dispersion only, super-gaussian shape, no loss, numerical

%we also compute the integration using fft, it is a numerical resolution
ws4 = 300;

%temporal domain operations
ts4 = linspace(-ws4*pw,ws4*pw,ns); % time values on which we compute the gaussian
% Eis4 = exp(-ts4.^2/(2*T0^2)); %gaussian values
Eis4 = exp(-0.5*(ts3/T0).^6); %super-gaussian input values
Ezs4 = zeros(1,ns); %pre allocation of the output field

%frequency domain operation

fs4 = ns/(2*ws4*pw); % fsample : sample at (number of indices/total window duration)
delta4 = 1/(2*ws4*pw); %interval between fft frequencies
fshift4 = (-ns/2:ns/2-1)*(fs4/ns); % zero-centered frequency range

% To use the fft function to convert the signal to the frequency domain,
% first identify a new input length that is the next power of 2 from the 
% original signal length. This will pad the signal X with trailing zeros 
% in order to improve the performance of fft.
% np = 2^nextpow2(ns); 
% 
% fs = np/(2*ws*pw); % fsample : sample at (number of indices/total window duration)
% delta = 1/(2*ws*pw); %interval between fft frequencies
% fshift = (-np/2:np/2-1)*(fs/np); % zero-centered frequency range

Fi4 = fftshift(fft(Eis4,ns));
Fi4 = abs(Fi4/ns);

% figure;
subplot(3,2,5);
plot(fshift4,Fi4);
title('Fourier Transform of the Super Gaussian to the fiber for third order');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

%integral using fft

% V1 = exp((1i/2)*beta2*length*fshift2.^2).*exp(-1i*fshift2*c);
for i = 1:ns
   Ezs4(i) = (1/(2*pi))*2*pi*delta4*sum(Fi4.*exp((1i/6)*beta3*length*fshift4.^3).*exp(-1i*fshift4*ts4(i)));
end
Ezs4 =  Ezs4/(delta4*(2*pi));

% figure;
subplot(3,2,6);
plot(ts4,abs(Ezs4));
title('Output signal out of the fiber for third order for Super Gaussian');
xlabel('time [s]');
ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');

