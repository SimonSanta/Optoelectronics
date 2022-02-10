%% Matlab script in order to simulate the laser diode.

%%
% 
% Instantion of all the constants.
clear all;
close all;

c = 299752458; % Speed of light
tau_p = 1*10^-12; % Photon life time
tau_r = 3*10^-9; % Radiative recombinaison time
tau_nr = 1*10^-9; % Nonradiative recombinaison time
tau_s = (1/tau_r + 1/tau_nr)^-1; % Carrier life time
beta = 1*10^-5; % Rate of spontaneous emission coupled to the laser mode
a = 3.0*10^-6; % Gain coefficient
n0 = 1*10^18; % Carrier density threshold
V = 1*10^-9; % Volume of cavity
e = 1.602*10^-19; % Elementary charge
lambda0 = 1.3*10^-6; % Usual operating wavelength of optical devices
f0 = c/lambda0; % Usual operating frequency of optical devices
alpha_int = 0.010; % Internal loss in m^-1 (10 cm^-1)
n = 3.5; % Reflection index of the active layer
eta = 1; % Quantum efficiency (100%)
l = 500*10^-6; % Cavity length
vg = c/n; %group celerity
hbar = (6.626*10^-34)/(2*pi); %normalised planck constant

%% First part : the static characteristics and parameters variations

% Instantion/preallocation of all the variable vectors.
A = zeros(1,10); % gain coefficient vector
R = zeros(1,10); % reflectivity vector
L = zeros(1,10); % cavity length vector

I = zeros(1,201); % current vector

Alpha_mirR = zeros(1,10); %loss vectors
Alpha_mirL = zeros(1,10);

Tau_pR =  zeros(1,10); % Photon life time vectors
Tau_pR =  zeros(1,10);

Pe0 = zeros(10,201);% power vector
PeA = zeros(10,201);% power vector
PeR = zeros(10,201);% power vector
PeL = zeros(10,201);% power vector


% At first we compute the value of power when the parameters?are not
% varying

Ith0 = (e*V/tau_s)*(n0+(1/(a*tau_p))); %threshold current with fixed parameters

r0 =((n-1)/(n+1))^2; % mirror reflectivity with n = 3.5
r1 = r0;
r2 = r0;
alpha_mir0 = (1/(2*l))*log(1/(r1*r2)); % L is the cavity length and internal loss in cm^-1

I = linspace(0,400,201)*10^-3; % current
Pe0 = ((hbar*2*pi*f0)/(2*e))*((1*alpha_mir0)/(alpha_mir0+alpha_int))*(I-Ith0);

% Then we compute the value of power when the parameters?are
% varying


R = linspace(r0*0.5,r0*1.5,10);
A = linspace(a*0.5,a*1.5,10);
L = linspace(l*0.5,l*1.5,10);


Alpha_mirR = (1/(2*l))*log(1./(R.*R)); %here R1 = R = R2
Alpha_mirL = (1./(2*L))*log(1./(r1*r2));

Tau_pR = (1./(vg*(Alpha_mirR+alpha_int)));
Tau_pL = (1./(vg*(Alpha_mirL+alpha_int)));

%the threshold current in each case
IthA = (e*V/tau_s)*(n0+(1./(A*tau_p))); %threshold current vector with varying a
IthR = (e*V/tau_s)*(n0+(1./(a*Tau_pR))); %threshold current vector with varying?R
IthL = (e*V/tau_s)*(n0+(1./(a*Tau_pL))); %threshold current vector with varying L

%the power profile in each case
PeA = ((hbar*2*pi*f0)/(2*e))*((1*alpha_mir0)./(alpha_mir0+alpha_int)).*(I-IthA(:));
PeR = ((hbar*2*pi*f0)/(2*e))*((1*Alpha_mirR(:))./(Alpha_mirR(:)+alpha_int)).*(I-IthR(:));
PeL = ((hbar*2*pi*f0)/(2*e))*((1*Alpha_mirR(:))./(Alpha_mirR(:)+alpha_int)).*(I-IthL(:));

%Let's set to zero all the values of the power which are negative
Pe0(Pe0<0)=0;
PeA(PeA<0)=0;
PeR(PeR<0)=0;
PeL(PeL<0)=0;

%plot of the graphs
figure;
subplot(2,2,1);
plot(I,Pe0(:,:));
title("P-I curve with fixed parameters");
xlabel("Current [A]");
ylabel("Power/Facet [W]");

subplot(2,2,2);
plot(I,PeA(:,:));
title("P-I curves with varying Gain coefficient a");
xlabel("Current [A]");
ylabel("Power/Facet [W]");
% p1 = plot(I,PeA(1,:));
% p10 = plot(I,PeA(10,:));
% legend([p1 p10],{'0.5*a','1.5*a'})

subplot(2,2,3);
plot(I,PeR(:,:));
xlim([0.2, 0.3]);
title("P-I curves with varying Reflectivity R");
xlabel("Current [A]");
ylabel("Power/Facet [W]");

subplot(2,2,4);
plot(I,PeL(:,:));
xlim([0.2, 0.3]);
title("P-I curves with varying Cavity length L");
xlabel("Current [A]");
ylabel("Power/Facet [W]");

%Variation of the threshold current

figure;
subplot(3,1,1);
plot(A,IthA);
title("Threshold current as a function of the gain a");
xlabel("gain a [cm^3/s]");
ylabel("Current [A]");

subplot(3,1,2);
plot(R,IthR);
title("Threshold current as a function of the Reflectivity R");
xlabel("Reflectivity R [.]");
ylabel("Current [A]");

subplot(3,1,3);
plot(L,IthL);
% xlim([0.2, 0.3]);
title("Threshold current as a function of the length L");
xlabel("length L [m]");
ylabel("Current [A]");

% %Variation of the ratio P/I
% figure;
% subplot(3,1,1);
% plot(A,IthA);
% title("Threshold current as a function of the gain a");
% xlabel("gain a [cm^3/s]");
% ylabel("Current [A]");
% 
% subplot(3,1,2);
% plot(R,IthR);
% title("Threshold current as a function of the Reflectivity R");
% xlabel("Reflectivity R [.]");
% ylabel("Current [A]");
% 
% subplot(3,1,3);
% plot(L,IthL);
% % xlim([0.2, 0.3]);
% title("Threshold current as a function of the length L");
% xlabel("length L [m]");
% ylabel("Current [A]");


%% Second part :Dynamic characteristics simulations
fwave = 100*10^6;
T = 1/fwave;

ti = linspace(0,2*T,300);
Ib = Ith0;
i = 2*Ib + Ib*square(2*pi*fwave*ti);
%i = 2*Ib;

% initial conditions
P0 = tau_p/(e*V)*(Ib-Ith0);
N0 = ((Ib/(e*V))+ a*n0*P0)/(1/tau_s + a*P0);

tspan = [0 2*T];
ic = [P0 N0];
options = odeset('RelTol',1e-8,'AbsTol',1e-8,'OutputFcn',@odeplot);

figure;
[t,f] = ode(@(t,f) myode(t,f,ti,i,a,beta,n0,tau_r,tau_p,e,V,tau_s), tspan, ic, options);

% Below f(1)=P and f(2)=N
function dfdt = myode(t,f,ti,i,a,beta,n0,tau_r,tau_p,e,V, tau_s)
h = interp1(ti,i,t);
dfdt = zeros(2,1);
dfdt(1) = a.*(f(2)-n0).*f(1) + (beta/tau_r)*f(2) - (1/tau_p)*f(1);
dfdt(2) = (1/e*V)*h - (1/tau_s)*f(2) - a.*(f(2)-n0).*f(1);
end


% function dydt = myode(t,y,ft,f,gt,g)
% f = interp1(ft,f,t); % Interpolate the data set (ft,f) at time t
% g = interp1(gt,g,t); % Interpolate the data set (gt,g) at time t
% dydt = -f.*y + g; % Evaluate ODE at time t



