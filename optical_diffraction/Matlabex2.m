
%% 1) Matlab script in order to simulate the faunhofer approximation.
% This code supports different aperture shape. The incident wave is a plane
% wave propagating in z direction, therefore the amplitude is constant over
% the plane (x,y) defined by the aperture.

%% Circular aperture case
%
% clear all;
close all;

lambda = 0.5*10^-6; % wavelength
apsize = 50*10^-6; % aperture size (10-100 micrometers)
% Faunhofer approximation is valid for dz >> apsize^2/lambda
dz = 1; % distance dz = 200*apsize^2/lambda
Ts = 1*lambda; %sampling period on the aperture

N = 512; % 512? 64? number of points in the aperture plane, better to be a power of 2
%N = 64;

%We sample the spatial domain
x = 1:N; %number of sample in x direction
y = 1:N; %number of sample in y direction
[X,Y] = meshgrid(x,y); %2 matrix of the x and y coordinates on each point of the 2D plane (x,y)
Xc = X-N/2-1;
Yc = Y-N/2-1;
dist = ((Xc.^2 + Yc.^2).^0.5)*Ts;

%We create a matrix with the values of the incident field at each point

E1 = zeros(N,N); %initialisation of the matrix

%?The values taken by this matrix depend on the shape of the aperture,
% in other word it depends on the domains on which the aperture is defined.
% If the we are inside the aperture's domain the value is 1 outside the
% value is zeros.


for i = 1:N
    for j = 1:N
        if dist(i,j) <= apsize/2 % if the point are inside de domain
            E1(i,j) = 1;
        end
    end
end

% We perform fft computations in 2D
F1 = zeros(N,N);
F1 = fftshift(fft2(fftshift(E1))); % the 2D dtft
F1 = abs(F1); % the amplitude of the values are taken
F1 = F1/max(max(F1)); % normalized with the max values
I1 = F1.^2;

SX = (X/N)*lambda*dz/Ts;
SY = (Y/N)*lambda*dz/Ts;

% Plots of the field F and the intensity pattern I

figure;
a1 = imagesc(E1);
colorbar;
saveas(a1,'ap1.png');

figure;
s11 = mesh(SX,SY,F1);
set(s11, 'EdgeColor', 'none');
s11.FaceColor = 'interp';
title('Normalized E in the observation plane for the circular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s11,'s11.png');

figure;
s12 = pcolor(SX,SY,F1);
set(s12, 'EdgeColor', 'none');
s12.FaceColor = 'interp';
title('Normalized E in the observation plane for the circular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
xlim([0.4 0.6]);
ylim([0.4 0.6]);
colorbar;
saveas(s12,'s12.png');

figure;
mesh(SX,SY,I1);
s13 = mesh(SX,SY,I1);
set(s13, 'EdgeColor', 'none');
s13.FaceColor = 'interp';
title('Normalized I pattern in the observation plane for the circular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s13,'s13.png');

figure;
s14 = pcolor(SX,SY,I1);
set(s14, 'EdgeColor', 'none');
s14.FaceColor = 'interp';
title('Normalized I pattern in the observation plane for the circular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
xlim([0.4 0.6]);
ylim([0.4 0.6]);
colorbar;
saveas(s14,'s14.png');

% imagesc(I1); title('Resulting intensity pattern distribution in the
% observation plane for the circular aperture'); xlabel('x-axis [m]');
% ylabel('y-axis [m]');

%% Rectangular aperture case

E2 = zeros(N,N); %initialisation of the matrix

%The values taken by this matrix depend on the shape of the aperture,
% in other word it depends on the domains on which the aperture is defined.
% If the we are inside the aperture's domain the value is 1 outside the
% value is zeros.

%Here, below, we define a new rectangular geometry for the aperture. This
%aperture's dimensions are the same as the original one in the x direction
%and in the y.

% lambda = 0.5*10^-6; apsize = 50*10^-6;

apsizex = apsize;
apsizey = apsize;

distx = Xc(1,:)*Ts; % here the indexes in x are transformed in spatial coordinates
disty = Yc(:,1)*Ts; % here the indexes in y are transformed in spatial coordinates

%The loop is modified to put the 1 insides the new aperture, 0 outside
for i = 1:N
    for j = 1:N
        if (((distx(i) <= apsizex/2) && (distx(i) >= -apsizex/2))&&((disty(j) <= apsizey/2) && (disty(j) >= -apsizey/2))) % if the point are inside de domain
            E2(i,j) = 1;
        end
    end
end

% We perform fft computations in 2D
F2 = zeros(N,N);
F2 = fftshift(fft2(fftshift(E2))); % the 2D dtft
F2 = abs(F2); % the amplitude of the values are taken
F2 = F2/max(max(F2)); % normalized with the max values
I2 = F2.^2;

% Plots of the field F and the intensity pattern I
figure;
a2 = imagesc(E2);
colorbar;
saveas(a2,'ap2.png');

figure;
s21 = mesh(SX,SY,F2);
set(s21, 'EdgeColor', 'none');
s21.FaceColor = 'interp';
title('Normalized E in the observation plane for the rectangular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s21,'s21.png');

figure;
s22 = pcolor(SX,SY,F2);
set(s22, 'EdgeColor', 'none');
s22.FaceColor = 'interp';
title('Normalized E in the observation plane for the rectangular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
xlim([0.4 0.6]);
ylim([0.4 0.6]);
colorbar;
saveas(s22,'s22.png');

% imagesc(F2); title('Resulting intensity pattern in the observation plane
% for the rectangular aperture'); xlabel('x-axis [m]'); ylabel('y-axis
% [m]');


figure;
s23 = mesh(SX,SY,I2);
set(s23, 'EdgeColor', 'none');
s23.FaceColor = 'interp';
title('Normalized I pattern in the observation plane for the rectangular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s23,'s23.png');

figure;
s24 = pcolor(SX,SY,I2);
set(s24, 'EdgeColor', 'none');
s24.FaceColor = 'interp';
title('Normalized I pattern in the observation plane for the rectangular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
xlim([0.4 0.6]);
ylim([0.4 0.6]);
colorbar;
saveas(s24,'s24.png');

% imagesc(I2); title('Resulting intensity pattern distribution in the
% observation plane for the rectangular aperture'); xlabel('x-axis [m]');
% ylabel('y-axis [m]');

%% Double slit case

E3 = zeros(N,N); %initialisation of the matrix

% The values taken by this matrix depend on the shape of the aperture, in
% other word it depends on the domains on which the aperture is defined. If
% the we are inside the aperture's domain the value is 1 outside the value
% is zeros.

%Here, below, we define a new double slit geometry for the aperture. This
%aperture's geometry is symetrical and parallel with respect to the y axis.



% lambda = 0.5*10^-6; apsize = 50*10^-6;

oneslitx = apsize/4;
gap = apsize/8;
apslity = apsize;
apslitx = 2*oneslitx + 2*gap;

distslitx = Xc(1,:)*Ts; % here the indexes in x are transformed in spatial coordinates
distslity = Yc(:,1)*Ts; % here the indexes in y are transformed in spatial coordinates

%The loop is modified to put the 1 insides the new aperture, 0 outside

for i = 1:N
    for j = 1:N
        if (((distslitx(i) <= apslitx/2) && (distslitx(i) >= -apslitx/2))&&((distslity(j) <= apslity/2) && (distslity(j) >= -apslity/2))) % if the point are inside de domain
            E3(i,j) = 1;
        end
    end
end

for i = 1:N
    for j = 1:N
        if (((distslitx(i) <= gap) && (distslitx(i) >= -gap))&&((distslity(j) <= apslity/2) && (distslity(j) >= -apslity/2))) % if the point are inside de domain
            E3(i,j) = 0;
        end
    end
end

% We perform fft computations in 2D
F3 = zeros(N,N);
F3 = fftshift(fft2(fftshift(E3))); % the 2D dtft
F3 = abs(F3); % the amplitude of the values are taken
F3 = F3/max(max(F3)); % normalized with the max values
I3 = F3.^2;

% Plots of the field F and the intensity pattern I

figure;
a3 = imagesc(E3);
colorbar;
saveas(a3,'ap3.png');

figure;
s31 = mesh(SX,SY,F3);
set(s31, 'EdgeColor', 'none');
s31.FaceColor = 'interp';
title('Normalized E in the observation plane for DS');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s31,'s31.png');

figure;
s32 = pcolor(SX,SY,F3);
set(s32, 'EdgeColor', 'none');
s32.FaceColor = 'interp';
title('Normalized E in the observation plane for DS');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
xlim([0.4 0.6]);
ylim([0.4 0.6]);
colorbar;
saveas(s32,'s32.png');

% imagesc(F3); title('Intensity pattern in the observation plane for DS');
% xlabel('x-axis [m]'); ylabel('y-axis [m]');


figure;
s33 = mesh(SX,SY,I3);
set(s33, 'EdgeColor', 'none');
s33.FaceColor = 'interp';
title('Normalized I pattern in the observation plane for DS');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s33,'s33.png');

figure;
s34 = pcolor(SX,SY,I3);
set(s34, 'EdgeColor', 'none');
s34.FaceColor = 'interp';
title('Normalized I pattern in the observation plane for for DS');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
xlim([0.4 0.6]);
ylim([0.4 0.6]);
colorbar;
saveas(s34,'s34.png');

% imagesc(I3); title('Intensity pattern in the observation plane for DS');
% xlabel('x-axis [m]'); ylabel('y-axis [m]');

%%  Periodic rectangular aperture case

E4 = zeros(N,N); %initialisation of the matrix

%The values taken by this matrix depend on the shape of the aperture,
% in other word it depends on the domains on which the aperture is defined.
% If the we are inside the aperture's domain the value is 1 outside the
% value is zeros.

%Here, below, we define a new rectangular geometry for the aperture. This
%aperture's dimensions are the same as the original one in the x direction
%and in the y.

% lambda = 0.5*10^-6; apsize = 50*10^-6;

apsizex = apsize/5;
apsizey = apsize/5;

distx = Xc(1,:)*Ts; % here the indexes in x are transformed in spatial coordinates
disty = Yc(:,1)*Ts; % here the indexes in y are transformed in spatial coordinates

inc =  max(distx)/10;

distx = distx - max(distx);
disty = disty - max(disty);

%The loop is modified to put the 1 insides the new aperture, 0 outside
%

for k = 1:19
    for g = 1:19
        distx = distx - max(distx) + inc*k;
        disty = disty - max(disty) + inc*g;
        for i = 1:N
            for j = 1:N
                if (((distx(i) <= apsizex/2) && (distx(i) >= -apsizex/2))&&((disty(j) <= apsizey/2) && (disty(j) >= -apsizey/2))) % if the point are inside de domain
                    E4(i,j) = 1;
                end
            end
        end
        
    end
end

% We perform fft computations in 2D
F4 = zeros(N,N);
F4 = fftshift(fft2(fftshift(E4))); % the 2D dtft
F4 = abs(F4); % the amplitude of the values are taken
F4 = F4/max(max(F4)); % normalized with the max values
I4 = F4.^2;

% Plots of the field F and the intensity pattern I

figure;
a4 = imagesc(E4);
colorbar;
saveas(a4,'ap4.png');

figure;
s41 = mesh(SX,SY,F4);
set(s41, 'EdgeColor', 'none');
s41.FaceColor = 'interp';
title('Normalized E in the observation plane for the periodic aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s41,'s41.png');

figure;
s42 = pcolor(SX,SY,F4);
set(s42, 'EdgeColor', 'none');
s42.FaceColor = 'interp';
title('Normalized E in the observation plane for the periodic aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
xlim([0.4 0.6]);
ylim([0.4 0.6]);
colorbar;
saveas(s42,'s42.png');

% imagesc(F4); title('Resulting intensity pattern in the observation plane
% for the periodic aperture'); xlabel('x-axis [m]'); ylabel('y-axis [m]');

% figure; s1 = pcolor(SX,SY,F4); set(s1, 'EdgeColor', 'none'); s1.FaceColor
% = 'interp';

figure;
s43 = mesh(SX,SY,I4);
set(s43, 'EdgeColor', 'none');
s43.FaceColor = 'interp';
title('Normalized I pattern in the observation plane for the periodic aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s43,'s43.png');

figure;
s44 = pcolor(SX,SY,I4);
set(s44, 'EdgeColor', 'none');
s44.FaceColor = 'interp';
title('Normalized I pattern in the observation plane for the periodic aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
xlim([0.4 0.6]);
ylim([0.4 0.6]);
colorbar;
saveas(s44,'s44.png');

% imagesc(I4); title('Resulting intensity pattern distribution in the
% observation plane for the periodic aperture'); xlabel('x-axis [m]');
% ylabel('y-axis [m]');

% figure; s2 = pcolor(SX,SY,I4); set(s2, 'EdgeColor', 'none'); s2.FaceColor
% = 'interp';

%% 2) Matlab script in order to simulate the Fresnel approximation.
% This code supports different aperture shape. The incident wave is a plane
% wave propagating in z direction, therefore the amplitude is constant over
% the plane (x,y) defined by the aperture.

%%  Circular aperture using Fresnel

% lambda = 0.5*10^-6; % wavelength apsize = 50*10^-6; % aperture size
% (10-100 micrometers) % Faunhofer approximation is valid for dz >>
% apsize^2/lambda dz = 1; % distance dz = 200*apsize^2/lambda Ts =
% 1*lambda; %sampling period on the aperture

%We sample the spatial domain of the observation plane

Rc = 70*lambda;
length = 300;
Circ = zeros (length,length);

%domain in observation plane?in polar coordinates
angle = linspace(0,2*pi,length); %number of sample for the angle
pau = linspace(0,apsize/2,length); %number of sample for the radius

%cartesian coordinates corresponding (one to one) to the polar ones for
%which we compute the field
dx = zeros(length,length);
dy = zeros(length,length);

for i = 1:1:length
    Circ(i,:) = circularshape(apsize,lambda,Rc, pau(i), angle(:));
    dx(i,:) = pau(i)*cos(angle(:));
    dy(i,:) = pau(i)*sin(angle(:));
end

Circ = Circ/max(max(Circ));
ICirc = Circ.^2;
% [Apau,Aangle] = meshgrid(pau,angle); ax = pau.*cos(theta); ay =
% pau.*sin(theta); [ax, ay] = pol2cart(angle,pau); [AX,AY] =
% meshgrid(ax,ay); dx = (dx/length)*lambda*R/Ts; dy =
% (dy/length)*lambda*R/Ts;

figure;
% surf(dx,dy,Circ);
s51 = mesh(dx,dy,Circ);
set(s51, 'EdgeColor', 'none');
s51.FaceColor = 'interp';
title('Normalized E in the observation plane using Fresnel and the circular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s51,'s51.png');

figure;
% s42 = pcolor(SX,SY,I4); set(s42, 'EdgeColor', 'none'); s42.FaceColor =
% 'interp'; title('Normalized I pattern in the observation plane for the
% periodic aperture'); xlabel('x-axis [m]'); ylabel('y-axis [m]');
% xlim([0.4 0.6]); ylim([0.4 0.6]);
s52 = pcolor(dx,dy,Circ);
set(s52, 'EdgeColor', 'none');
s52.FaceColor = 'interp';
title('Normalized E in the observation plane using Fresnel and the circular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s52,'s52.png');

% xlim([0.4 0.6]); ylim([0.4 0.6]); surfc(dx,dy,Circ); hold on
% imagesc(dx,dy,Circ);

figure;
% surf(dx,dy,ICirc);
s53 = mesh(dx,dy,ICirc);
set(s53, 'EdgeColor', 'none');
s53.FaceColor = 'interp';
title('Normalized I in the observation plane using Fresnel and the circular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s53,'s53.png');

figure;
s54 = pcolor(dx,dy,ICirc);
set(s54, 'EdgeColor', 'none');
s54.FaceColor = 'interp';
title('Normalized I in the observation plane using Fresnel and the circular aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s54,'s54.png');

% surfc(dx,dy,ICirc); hold on imagesc(dx,dy,ICirc);

%%  knife edge aperture using Fresnel

Rkf = 800*lambda;
length2 = 100;
KF = zeros (2*length2+1,2*length2+1);

%we display the aperture shape
disc = 200; % number of element in each dimesion
S = zeros(disc+1,disc+1); % matrix of the aperture

% xkf = linspace(0,apsize/2,disc/2+1);
xkf = linspace(-2*apsize/2,2*apsize/2,disc+1);
ykf = linspace(-2*apsize/2,2*apsize/2,disc+1); % discretisation of the aperture domain
[Sx,Sy] = meshgrid(xkf,ykf);

for gy = 1:1:disc+1 % iterate on the radius
    for hx = 1:1:disc+1 % iterate on the angle
        if(xkf(hx) <= 0)
            S(gy,hx)=0;
        else
            S(gy,hx)=1;
        end
    end
end

figure;
Sshape = pcolor(Sx,Sy,S);
set(Sshape, 'EdgeColor', 'none');
Sshape.FaceColor = 'interp';
title('KF aperture shape');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(Sshape,'Sshape.png');

%cartesian coordinates corresponding to the point in 2D in which we compute
%the resulting field, it is the observation plane

xprime = linspace(-apsize/4,2*apsize/2,2*length2+1);
yprime = linspace(-2*apsize/2,2*apsize/2,2*length2+1);

for i = 1:1:2*length2+1
    KF(i,:) = kfshape(apsize, lambda,Rkf, xprime(:), yprime(i));
    i
end

KF = KF./max(max(KF));
IKF = KF.^2;

[Xp,Yp] = meshgrid(xprime,yprime);

% graphes still to do
figure;
% surf(dx,dy,Circ);
s61 = mesh(Xp,Yp,KF);
set(s61, 'EdgeColor', 'none');
s61.FaceColor = 'interp';
title('Normalized E in the observation plane using Fresnel and the KF aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s61,'s61.png');

figure;
% s42 = pcolor(SX,SY,I4); set(s42, 'EdgeColor', 'none'); s42.FaceColor =
% 'interp'; title('Normalized I pattern in the observation plane for the
% periodic aperture'); xlabel('x-axis [m]'); ylabel('y-axis [m]');
% xlim([0.4 0.6]); ylim([0.4 0.6]);
s62 = pcolor(Xp,Yp,KF);
set(s62, 'EdgeColor', 'none');
s62.FaceColor = 'interp';
title('Normalized E in the observation plane using Fresnel and the KF aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s62,'s62.png');

% xlim([0.4 0.6]); ylim([0.4 0.6]); surfc(dx,dy,Circ); hold on
% imagesc(dx,dy,Circ);

figure;
% surf(dx,dy,ICirc);
s63 = mesh(Xp,Yp,IKF);
set(s63, 'EdgeColor', 'none');
s63.FaceColor = 'interp';
title('Normalized I in the observation plane using Fresnel and the KF aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s63,'s63.png');

figure;
s64 = pcolor(Xp,Yp,IKF);
set(s64, 'EdgeColor', 'none');
s64.FaceColor = 'interp';
title('Normalized I in the observation plane using Fresnel and the KF aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(s64,'s64.png');


%% DS aperture using Fresnel

%
Rds = 100*lambda;
length3 = 100;
DS = zeros (length3+1,length3+1);

% Let's draw the aperture domain
% vectors for the aperture domain
sample = 200;

%space between the two parallel slits
gapx = apsize/10;
gapy = 0;

%dimensions of the slits
slitx = apsize/10;
slity = 3*apsize/10;

% total domain on which we draw the aperture
xds = linspace(-2*apsize/2,2*apsize/2,sample+1);
yds = linspace(-2*apsize/2,2*apsize/2,sample+1);
D = zeros(sample+1,sample+1);
[Dx,Dy] = meshgrid(xds,yds);

%apertures limits in the space
xmin = gapx/2;
xmax = gapx/2 + slitx;
ymin = gapy/2;
ymax = gapy/2 + slity;

% we fill the space where the aperture is located
for gy = 1:1:sample+1 % iterate on the radius
    for hx = 1:1:sample+1 % iterate on the angle
        if((yds(gy) >= -ymax && yds(gy) <= -ymin) || ...
                (yds(gy) <= ymax && yds(gy) >= ymin))
            
            if((xds(hx) >= -xmax && xds(hx) <= -xmin) || ...
                    (xds(hx) <= xmax && xds(hx) >= xmin))
                D(gy,hx) = 1;
            else
                D(gy,hx) = 0;
            end
            
        else
            D(gy,hx) = 0;
        end
    end
end

% Aperture drawing
figure;
Dshape = pcolor(Dx,Dy,D);
set(Dshape, 'EdgeColor', 'none');
Dshape.FaceColor = 'interp';
title('Normalized E in the observation plane using Fresnel and the DS aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
saveas(Dshape,'Dshape.png');

%cartesian coordinates corresponding to the point in 2D in which we compute
%the resulting field, it is the observation plane

xdsbis = linspace(-apsize/2,apsize/2,length3+1);
ydsbis = linspace(-apsize/2,apsize/2,length3+1);
[Dxbis,Dybis] = meshgrid(xdsbis,ydsbis);

for k = 1:1:length3+1
    DS(k,:) = slitshape(apsize, lambda,Rds,xmax,xmin,ymax,ymin,xdsbis(:),ydsbis(k));
    k
end

DS = DS./max(max(DS));
IDS = DS.^2;

figure;
% surf(dx,dy,Circ);
s71 = mesh(Dxbis,Dybis,DS);
set(s71, 'EdgeColor', 'none');
s71.FaceColor = 'interp';
title('Normalized E in the observation plane using Fresnel and the DS aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
xlim([-2.5 2.5]*10^-5);
ylim([-2.5 2.5]*10^-5);
saveas(s71,'s71.png');

figure;
% s42 = pcolor(SX,SY,I4); set(s42, 'EdgeColor', 'none'); s42.FaceColor =
% 'interp'; title('Normalized I pattern in the observation plane for the
% periodic aperture'); xlabel('x-axis [m]'); ylabel('y-axis [m]');
% xlim([0.4 0.6]); ylim([0.4 0.6]);
s72 = pcolor(Dxbis,Dybis,DS);
set(s72, 'EdgeColor', 'none');
s72.FaceColor = 'interp';
title('Normalized E in the observation plane using Fresnel and the DS aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
xlim([-2.5 2.5]*10^-5);
ylim([-2.5 2.5]*10^-5);
saveas(s72,'s72.png');
% surfc(dx,dy,Circ); hold on
% imagesc(dx,dy,Circ);

figure;
% surf(dx,dy,ICirc);
s73 = mesh(Dxbis,Dybis,IDS);
set(s73, 'EdgeColor', 'none');
s73.FaceColor = 'interp';
title('Normalized I in the observation plane using Fresnel and the DS aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
xlim([-2.5 2.5]*10^-5);
ylim([-2.5 2.5]*10^-5);
colorbar;
saveas(s73,'s73.png');

figure;
s74 = pcolor(Dxbis,Dybis,IDS);
set(s74, 'EdgeColor', 'none');
s74.FaceColor = 'interp';
title('Normalized I in the observation plane using Fresnel and the DS aperture');
xlabel('x-axis [m]');
ylabel('y-axis [m]');
colorbar;
xlim([-2.5 2.5]*10^-5);
ylim([-2.5 2.5]*10^-5);
saveas(s74,'s74.png');

%% Part of the code dedicated to matlab function designed to solve the problem in the different cases

%% 1)
%Function used to compute the integral in the Fresnel's case for the
%circular aperture, this integral is calculated in terms of polar
%coordinates
function field = circularshape(apsize, lambda,z, pau, angle)

sum = 0;
k = 2*pi/lambda;

disc = 100; % the number of values for the discretization of the aperture
theta = linspace(0,2*pi,disc+1); % angle of the aperture domain
r = linspace(0,apsize/2,disc+1); % radius of the aperture domain

for g = 1:1:disc % iterate on the radius
    for h = 1:1:disc % iterate on the angle
        
        arg = sqrt(r(g)^2 + z^2 + pau^2 - 2*r(g)*pau*cos(theta(h)));
        sum = sum + (exp(-1i*k*arg)/arg);
        
    end
end

field = abs(sum);

end

%% 2)
%Function used to compute the integral in the Fresnel's case for the
%knife-edge aperture, this integral is calculated in terms cartesian
%coordinates

function field = kfshape(apsize, lambda,z, xprime, yprime)

sum = 0;
% k = 2*pi/(z*lambda);
k = 2*pi/(lambda);
% k= k/z;

disc = 200; % the number of values for the discretization of the aperture
% xkf = linspace(0,apsize/2,disc/2+1); ykf =
% linspace(-apsize/2,apsize/2,disc+1); % angle of the aperture domain

xkf = linspace(-2*apsize/2,2*apsize/2,disc+1); %to have good shapes replace by 4
ykf = linspace(-2*apsize/2,2*apsize/2,disc+1); % aperture domain 201 points in x ad y

for gy = 1:1:disc+1 % iterate on y
    for hx = 1:1:disc+1 % iterate on x
        if(xkf(hx) <= 0)
            sum = sum + 0;
        else
            arg = (xprime - xkf(hx)).^2 + (yprime - ykf(gy))^2;
            argd = arg./sqrt(z^2 + xprime.^2 + yprime^2 - xkf(hx)^2 - ykf(gy)^2);
            div = 1./sqrt(z^2 + xprime.^2 + yprime^2 - xkf(hx)^2 - ykf(gy)^2);
            %below the second method with the atan
            %             arg = sqrt(xkf(hx)^2 + ykf(gy)^2 + xprime.^2 + yprime^2 + z^2 ...
            %                     - 2*(xkf(hx)^2 + ykf(gy)^2)*(xprime.^2 + yprime^2).* ...
            %                     cos(atan(xprime./yprime) - atan(ykf(gy)/xkf(hx))));
            %
            
            %below the first method without the atan
            %             arg = sqrt(xkf(hx)^2 + ykf(gy)^2 + xprime.^2 + yprime^2 + z^2 ...
            %                     - 2*(xkf(hx)^2 + ykf(gy)^2)*(xprime.^2 + yprime^2).* ...
            %                     cos(0 - atan(ykf(gy)/xkf(hx))));
            
            %             sum = sum + (exp(-1i*k*arg)./arg);
            sum = sum + (div).*(exp(-1i*k*argd));
            %            sum = sum + (exp(-1i*k*argd));
        end
    end
end

field = abs(sum);

end

%% 3)
%Function used to compute the integral in the Fresnel's case for the
%DS aperture, this integral is calculated in terms cartesian
%coordinates
function field = slitshape(apsize,lambda,z,xmax,xmin,ymax,ymin,xprime,yprime)

sum = 0;
k = 2*pi/lambda;

disc = 100; % the number of values for the discretization of the aperture
% xkf = linspace(0,apsize/2,disc/2+1);
% ykf = linspace(-apsize/2,apsize/2,disc+1);

xds = linspace(-apsize/2,apsize/2,disc+1);
yds = linspace(-apsize/2,apsize/2,disc+1);

for gy = 1:1:disc+1 % iterate on the radius
    for hx = 1:1:disc+1 % iterate on the angle
        if((yds(gy) >= -ymax && yds(gy) <= -ymin) || ...
                (yds(gy) <= ymax && yds(gy) >= ymin))
            
            if((xds(hx) >= -xmax && xds(hx) <= -xmin) || ...
                    (xds(hx) <= xmax && xds(hx) >= xmin))
                
                arg = (xprime - xds(hx)).^2 + (yprime - yds(gy))^2;
                argd = arg./sqrt(z^2 + xprime.^2 + yprime^2 - xds(hx)^2 - yds(gy)^2);
                div = 1./sqrt(z^2 + xprime.^2 + yprime^2 - xds(hx)^2 - yds(gy)^2);
                
                %                 arg = sqrt(xds(hx)^2 + yds(gy)^2 + xprime.^2 + yprime^2 + z^2 ...
                %                 - 2*(xds(hx)^2 + yds(gy)^2)*(xprime.^2 + yprime^2).*...
                %                 cos(atan(xprime./yprime) - atan(yds(gy)/xds(hx))) );
                
                %                 arg = sqrt(xds(hx)^2 + yds(gy)^2 + xprime.^2 + yprime^2 + z^2 ...
                %                     - 2*(xds(hx)^2 + yds(gy)^2)*(xprime.^2 + yprime^2).*...
                %                     cos(0 - atan(yds(gy)/xds(hx))) );
                %                  sum = sum + (exp(-1i*k*arg)./arg);
                
                sum = sum + (div).*(exp(-1i*k*argd));
            else
                sum = sum + 0;
            end
            
        else
            sum = sum + 0;
        end
    end
end


field = abs(sum);

end
