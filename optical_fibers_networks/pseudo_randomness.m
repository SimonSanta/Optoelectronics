
close all;

% constants and parameters all units in SI units
D = (20*10^-12)/(10^-9*10^3); %dispersion parameter : s/m^2
S = (0.07*10^-12)/(10^-18*10^3); %dispersion slope : s/m^3
lambda0 = 1500*10^-9; %wavelength : m
c = 299792458; %speed of light : m/s
f0 = c/lambda0;
w0 = 2*pi*f0;
% pw = 50*10^-12; %pulse width : s
length = 100*1000; %length of fiber : m
beta2 = -D*lambda0^2/(2*pi*c); %propagation constant of second order
beta3 = S*(lambda0^2/(2*pi*c))^2; %propagation constant of third order
% pw = 0.1; %pulse width : s
pw = 1*10^-12; %pulse width : s
T0 = pw/1.665; % pulse width T0

ns = 1000;
ws = 100;

%% first step, the pseudo random bit sequence using n = 15
nrbs = 15;
yseq1 = ltePRBS(162,nrbs)';
yseq1 = double(yseq1);
xseq1 = linspace(-(nrbs/2),(nrbs/2),nrbs);

%plot1
figure;
stem(xseq1,yseq1);
xlabel('time [s]');
ylabel('Normalized amplitude');

%% second step, the gaussian filter
xseq2 = linspace(-(nrbs/2)*ws*pw,(nrbs/2)*ws*pw,nrbs*ns);
Gausshape = exp(-xseq2.^2/(2*T0^2));

% xgaussbit = linspace(-(nrbs/2)-1,(nrbs/2)+1,nrbs*100);

%plot2
figure;
plot(xseq2,Gausshape); % drawing of an isolated gaussian signal
xlabel('time [s]');
ylabel('Normalized amplitude');

%% Here we implement the code to put the delta at the places of the speudo random bits

% Below, vector setting if we have to place a gaussian shape or not

xbin = linspace(-(nrbs/2+1),(nrbs/2+1),nrbs+2);
ybin = ltePRBS(162,nrbs)'; %again the bin sequence?this time just to decide if we draw the signal or not
ybin2 = zeros(1,17);
ybin2(1,2:16) = ybin;

%plot3
figure;
stem(xbin,ybin2);
xlabel('time [s]');
ylabel('Normalized amplitude');

%%%%%%%%%%%%%%%
%We will, at the same time as we put the gaussian, compute the effect of
%the dispersion on those

%Each pulse is isolated and occupies one row of the ouputfield matrix

Ezgauss1 = zeros(nrbs+2,nrbs*ns); %pre allocation of the output field, matrix with nrbs rows of same size as the input field Egaussbit
Ezgauss2 = zeros(nrbs+2,nrbs*ns); %third order dispersion
% below parameters usefull for the dispersion effect

fs1 = (nrbs*ns)/(2*ws*pw); % fsample : sample at (number of indices/total window duration)
% delta = 1/(2*ws*pw); %interval between fft frequencies
deltagauss1 = 1/(2*ws*pw); %interval between fft frequencies
fshift1 = (-(nrbs*ns)/2:(nrbs*ns)/2-1)*(fs1/(nrbs*ns)); % zero-centered frequency range

% We switch to the frequencies
% Fgaussi = fftshift(fft(Egaussbit,nrbs*ns));
% % Fi = fft(Eis,ns)
% Fgaussi = abs(Fgaussi/nrbs*ns);

%%%%%%%%%%%%%%%

%For loop putting the gaussian
Egaussbit  = zeros(1,nrbs*ns);
Egaussbitprime  = zeros(1,nrbs*ns);
expx = zeros(1,nrbs*ns);
Fexp = zeros(1,nrbs*ns);
g = 1;

for i = 1:nrbs+2
    if(ybin2(i))
%         xseq2 = linspace(-(nrbs/2)*ws*pw,(nrbs/2)*ws*pw,nrbs*ns);
        xfor = xseq2 + (nrbs/2+1)*ws*pw - (i)*ws*pw*((nrbs+2)/(nrbs+1));
        Egaussbit = Egaussbit + exp(-xfor.^2/(2*T0^2));
        
        expx = exp(-xfor.^2/(2*T0^2));
        Fexp = fftshift(fft(expx,nrbs*ns));
        Fexp = abs(Fexp/(nrbs*ns));
        
        %%%%%%%%%%%%%
        % below, the operation simulating the effect of the second order dispersion
        for k = 1:nrbs*ns
            Ezgauss1(g,k) = (1/(2*pi))*2*pi*deltagauss1*sum(Fexp.*exp((1i/2)*beta2*length*fshift1.^2).*exp(-1i*fshift1*xfor(k)));
            Ezgauss2(g,k) = (1/(2*pi))*2*pi*deltagauss1*sum(Fexp.*exp((1i/6)*beta3*length*fshift1.^3).*exp(-1i*fshift1*xfor(k)));
        end
 
        Ezgauss1(g,:) =  Ezgauss1(g,:)./(deltagauss1*(2*pi));
        Ezgauss2(g,:) =  Ezgauss2(g,:)./(deltagauss1*(2*pi));
        %%%%%%%%%%%%%%
        g = g+1;
    else
        Ezgauss1(g,:) =  0;
        Ezgauss2(g,:) = 0;
        g = g+1;
    end
end

% Egaussbitprime = sum(Ezgauss1);
tgaussbit = linspace((-(nrbs/2)-1)*ws*pw,((nrbs/2)+1)*ws*pw,nrbs*ns);

%plot4
figure;
plot(tgaussbit,Egaussbit);
xlabel('time [s]');
ylabel('Normalized amplitude');

%plot5

Ezgauss1(2,7500:15000) = 0;
Ezgauss1(14,1:7500) = 0;
Ezgauss1(13,1:7500) = 0;
Ezgauss2(2,7500:15000) = 0;
Ezgauss2(14,1:7500) = 0;
Ezgauss2(13,1:7500) = 0;

figure;
plot(tgaussbit,abs(Ezgauss1(:,:)));
xlabel('time [s]');
ylabel('Normalized amplitude');

%plot6

figure;
plot(tgaussbit,abs(Ezgauss2(:,:)));
xlabel('time [s]');
ylabel('Normalized amplitude');
%% fft computation

% ws*pulsewidth represent the number

% below, the usefull parameters
% pw = 50*10^-12; %pulse width : s
% T0 = pw/1.665; % pulse width T0
%
% ns = 1000;
% ws = 100;

% below, the usefull quantities. The second is the input signal to the
% fiber, and represent a pseudo random bit sequence of gaussian pulses.

% tgaussbit
% Egaussbit
%% analitycal method

% %we compute the output field for second order dispersion here the Fourier
% %transform of the input was computed analytically
% Ez1 = zeros(1,nrbs*ns);
% int1 = zeros(1,nrbs*ns);
% k = 2*T0^2;
% 
% fun1 = @(dw,c) sqrt(k*pi)*exp(-pi^2*k*(dw.^2)).*exp((1i/2)*beta2*length*dw.^2).*exp(-1i*dw*c);
% for i = 1:nrbs*ns
%     int1(i) = integral(@(dw)fun1(dw,tgaussbit(i)),-inf,inf);
%     i
% end
% Ez1 = (1/(2*pi))*int1;
% 
% %plot6
% figure;
% plot(tgaussbit,Ez1);
% xlabel('time [s]');
% ylabel('Normalized amplitude of the pulse Ein/Einmax');


%% first step, the pseudo random bit sequence using n = 15

% Ezgauss = zeros(1,nrbs*ns); %pre allocation of the output field, same size as the input field Egaussbit
% 
% 
% fs = (nrbs*ns)/(2*ws*pw); % fsample : sample at (number of indices/total window duration)
% % delta = 1/(2*ws*pw); %interval between fft frequencies
% deltagauss = 1/(2*ws*pw); %interval between fft frequencies
% fshift = (-(nrbs*ns)/2:(nrbs*ns)/2-1)*(fs/(nrbs*ns)); % zero-centered frequency range
% 
% % We switch to the frequencies
% Fgaussi = fftshift(fft(Egaussbit,nrbs*ns));
% % Fi = fft(Eis,ns)
% Fgaussi = abs(Fgaussi/(nrbs*ns));
% 
% %plot7
% figure;
% plot(fshift,Fgaussi);
% title('Fourier Transform of the input signal to the fiber for second order- Gaussian - Numeric');
% xlabel('time [s]');
% ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');
% 
% % below, the operation simulating the effect of the second order dispersion
% for j = 1:nrbs*ns
%     Ezgauss(j) = (1/(2*pi))*2*pi*deltagauss*sum(Fgaussi.*exp((1i/2)*beta2*length*fshift.^2).*exp(-1i*fshift*tgaussbit(j)));
%     
% end
% Ezgauss =  Ezgauss/(deltagauss*(2*pi));
% 
% %plot8
% figure;
% plot(tgaussbit,abs(Ezgauss));
% title('Output signal out of the fiber for second order - Gaussian - Numerical');
% xlabel('time [s]');
% ylabel('Normalized amplitude of the pulse Eout/Einmax[V/m]');
% 
% 
% %third step, the convolution between them
% %
% %  yconv = conv(yseq1,Gausshape);
% %  convlength = nrbs*100+nrbs-1;
% %  xconv = linspace(-(nrbs),(nrbs),convlength);
% 
% %plot3
% % figure;
% % plot(xconv,yconv);
% 
% 
% % xseq2 = linspace(-(nrbs/2),(nrbs/2),nrbs*10);
% % yseq2 = interp1(xseq1,yseq1,xseq2,'gaussian');
% %
% % figure;
% % stem(xseq2,yseq2)