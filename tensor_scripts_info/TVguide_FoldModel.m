clear all
close all
%%
% Reads in data from xls file, with sample # in first column, then density in
% g/cm^3, then the 21 tensor components:
% Samp#,rho, C11, C12, C13, C14, C15, C16, C22, C23,...

Data = xlsread('PelonaDB_XYZorientation.xlsx');
%% Choosing which tensor to plot from the DB
% change this to choose your sample
SampleN = 70;
FigName = ['Tensor # ',num2str(SampleN)];


Samps = find(Data(:,1) == SampleN); 
Tensors = Data(Samps,3:23);
Rhos = Data(Samps,2);
Samples = num2str(SampleN);

C = Tensors;
C = [C(1:6);...
    C(2), C(7:11);...
    C(3), C(8), C(12:15);...
    C(4), C(9), C(13), C(16:18);...
    C(5), C(10), C(14), C(17), C(19:20);...
    C(6), C(11), C(15), C(18), C(20), C(21)];
rho = Rhos;

%% Fold model is currently made using a sine function. Fold axis is N-S, and 
% the sine function specifies the dip. We could put other shapes in if we
% wanted, but this works well, and is easy. The figures that come out are a
% cross section showing the shape of the fold, and stereonets of the averaged tensor that are
% geographic xyz, so fold axis is N-S.

% Angles of the sine function
Ang = (-180:180);

% Amplitude
Amp = 25;

% Period
Per = 4;

% Resolution, how many tensors to put in
ResN = 100;

% Model fold
Shape = Amp*sind(Per*Ang);

% Tensor locations
TLoc = round(linspace(1,length(Ang)-1,ResN));

% Dip angles taken from derivative of sine curve (or whatever Shape is
% defined)
dip = real(-1*atand(diff(Shape)));
dip = dip(TLoc);

% For a 2D fold cross-section the strike is always into the page, needs to be specified
% for the rotation, so setting it to zero
strike = zeros(length(dip),1);

% Trend is the direction of the lineation, 0 would be parallel to fold
% hinge, 90 is perpendicular to fold hinge
trend = 90*ones(length(dip),1);

%% Cross-section plot
figure(1)
plot(Ang,Shape)
hold on
plot(Ang(TLoc),Shape(TLoc),'kd')
axis equal
plot(Ang(TLoc),dip)
legend({'Structure';'Tensor location';'dip angle'})
%% Empty matrix for all rotated tensors to go into
RotatedCij = zeros(length(dip(:,1)),22);
%%
% If inputting strike, dip, and trend to specify a rotation, need to
% calculate rake.
rake = zeros(length(dip),1);
for i = 1:length(dip)
    if trend(i) - strike(i) == 180;
        rake(i) = 180;
    elseif trend(i) < strike(i)
        trend(i) = trend(i) + 360;
    end
    if dip == 0 % When dip is 0, rake is equal to trend
        rake(i) = trend(i);
    elseif dip == 90 % When dip is 90, trend can only be parallel to strike
        rake(i) = plunge(i);
    else
        SL = cosd(trend(i) - strike(i));
        L1 = sind(trend(i)-strike(i));
        L2 = L1/cosd(dip(i));
        rake(i) = atand(L2/SL);
        if rake(i) < 0
            rake(i) = 180 + rake(i);
        end
    end
    rake(i) = -1*rake(i);
end
%% These are the specific rotations around each axis

for i = 1:length(dip)
    ya(i) = dip(i); % Rotates around top, around Y if XY orientation, DIP
    xa(i) = 0; % Rotates around right, around X ''
    za(i) = strike(i); % Rotates around center, around Z '', STRIKE
    yrake(i) = rake(i)-90; % This also rotates around center, before the dip rotation, to put lineation at the right rake
end

%% Creating the rotation matrices and applying the rotation

for s = 1:length(dip)
        
        zr = -1*za(s)*pi/180;
        xr = -1*xa(s)*pi/180;
        yr = 1*ya(s)*pi/180;
        yrk = 1*yrake(s)*pi/180;

        Rx = [1 0 0 0 0 0;...
            0 cos(xr)^2 sin(xr)^2 2*cos(xr)*sin(xr) 0 0;...
            0 sin(xr)^2 cos(xr)^2 -2*cos(xr)*sin(xr) 0 0;...
            0 -cos(xr)*sin(xr) cos(xr)*sin(xr) cos(xr)^2-sin(xr)^2 0 0;...
            0 0 0 0 cos(xr) -sin(xr);...
            0 0 0 0 sin(xr) cos(xr)];
        Ry = [cos(yr)^2 0 sin(yr)^2 0 2*cos(yr)*sin(yr) 0;...
            0 1 0 0 0 0;...
            sin(yr)^2 0 cos(yr)^2 0 -2*cos(yr)*sin(yr) 0;...
            0 0 0 cos(yr) 0 -sin(yr);...
            -cos(yr)*sin(yr) 0 cos(yr)*sin(yr) 0 cos(yr)^2-sin(yr)^2 0;...
            0 0 0 sin(yr) 0 cos(yr)];
        Rz = [cos(zr)^2 sin(zr)^2 0 0 0 -2*cos(zr)*sin(zr);...
            sin(zr)^2 cos(zr)^2 0 0 0 2*cos(zr)*sin(zr);...
            0 0 1 0 0 0;...
            0 0 0 cos(zr) sin(zr) 0;...
            0 0 0 -sin(zr) cos(zr) 0;...
            cos(zr)*sin(zr) -cos(zr)*sin(zr) 0 0 0 cos(zr)^2-sin(zr)^2];
        Ryrake = [cos(yrk)^2 sin(yrk)^2 0 0 0 -2*cos(yrk)*sin(yrk);...
            sin(yrk)^2 cos(yrk)^2 0 0 0 2*cos(yrk)*sin(yrk);...
            0 0 1 0 0 0;...
            0 0 0 cos(yrk) sin(yrk) 0;...
            0 0 0 -sin(yrk) cos(yrk) 0;...
            cos(yrk)*sin(yrk) -cos(yrk)*sin(yrk) 0 0 0 cos(yrk)^2-sin(yrk)^2];
        
       % Apply rotations in order: 1. Rotate lineation to correct rake
       % anlge, 2. rotate to correct dip, 3. rotate to correct strike
        CR = C;
        CR = Ryrake*CR*Ryrake'; % Rotates to rake
        CR = Rx*CR*Rx'; % This doesn't actually do anything because no rotation angle around x
        CR = Ry*CR*Ry'; % Rotates to dip
        CR = Rz*CR*Rz'; % Rotates to strike
        
        
        
        RotatedCij(s,1) = rho;
        RotatedCij(s,2:7) = CR(1,:);
        RotatedCij(s,8:12) = CR(2,2:6);
        RotatedCij(s,13:16) = CR(3,3:6);
        RotatedCij(s,17:19) = CR(4,4:6);
        RotatedCij(s,20:21) = CR(5,5:6);
        RotatedCij(s,22) = CR(6,6);

end
%% Averaging V is for Voigt (using Cij), R is for Reuss (using inverse of Cij)
TensorsV = cell(length(dip),1);
TensorsR = cell(length(dip),1);
V_ave = zeros(6,6);
R_ave = zeros(6,6);

for i = 1:ResN
    C = RotatedCij(i,2:end);
    TensorsV{i} = [C(1:6);...
        C(2), C(7:11);...
        C(3), C(8), C(12:15);...
        C(4), C(9), C(13), C(16:18);...
        C(5), C(10), C(14), C(17), C(19:20);...
        C(6), C(11), C(15), C(18), C(20), C(21)];
    
    TensorsR{i} = inv(TensorsV{i});
end

% This loop averages the tensor components for the V and R averages
for i = 1:36
    PV = zeros(length(dip),1);
    PR = zeros(length(dip),1);
    for j = 1:ResN
        PV(j) = TensorsV{j}(i);
        PR(j) = TensorsR{j}(i);
    end
    V_ave(i) = mean(PV);
    R_ave(i) = mean(PR);
end
% VRH average is the average of V and R stiffness tensors
VRH = (V_ave + inv(R_ave))/2;
%% This is just the Christoffel and plotting codes, which wouldn't need to be repeated here necessarily

Res = 3;
for s = 1
C = VRH;

rho = rho;
%%
%   mapping 4th-order Cijkl indices to Voigt notation 
Voigt_ijkl = [ 1, 6, 5 ; 6, 2, 4 ; 5, 4, 3 ] ; 

[Xx,Yy,Zz] = sphere(360/Res);

for phi = 1:size(Xx,1)

     phi_index = phi;

    for theta = 1:size(Xx,1)
         theta_index = theta;
        T = zeros(3, 3) ; 
        TT = zeros(3, 3) ; 
        
        %   direction cosines of this unit vector
        X = [ Xx(theta,phi), Yy(theta,phi), Zz(theta,phi) ] ;
        
        %   calculate Tik, Christoffel tensor 
        for i = 1:3 
            for k = 1:3 
                T(i, k) = 0.0 ; 
                for j = 1:3 
                    for l = 1:3 
                        m = Voigt_ijkl(i, j) ; 
                        n = Voigt_ijkl(k, l) ;
                        T(i, k) = T(i, k) + C(m, n) * X(j) * X(l) ; 
                    end ; 
                end ; 
            end ; 
        end ; 
        
        %   form TTij = Tik * Tjk
        for i = 1:3 
            for j = 1:3 
                TT(i, j) = 0.0 ; 
                for k = 1:3 
                    TT(i, j) = TT(i, j) + T(i, k) * T(j, k) ; 
                end ; 
            end ; 
        end ; 
        
        %   get eigenvalues of TTij 
        [eigenVecs, eigenVals] = eig(TT) ; 
        velp(theta, phi) = sqrt( sqrt( eigenVals(3,3) ) / rho )*10 ; %Vp as a function of theta & phi
        vels1(theta, phi) = sqrt( sqrt( eigenVals(2,2) ) / rho )*10 ; %Vs1 as a function of theta & phi
        vels2(theta, phi) = sqrt( sqrt( eigenVals(1,1) ) / rho )*10 ; %Vs2 as a function of theta & phi
    
        VpPolxyz(theta,phi,:) = eigenVecs(:,3)';
        Vs1Polxyz(theta,phi,:) = eigenVecs(:,2)';
        Vs2Polxyz(theta,phi,:) = eigenVecs(:,1)';
       
    end ; 
    
end ; 

%% Isotropic properties
Cij = C;
Sij = inv(Cij);

Kvoigt = (1/9)*(Cij(1,1)+Cij(2,2)+Cij(3,3)+2*(Cij(1,2)+Cij(1,3)+Cij(2,3)));
Kreuss = 1/(Sij(1,1)+Sij(2,2)+Sij(3,3)+2*(Sij(1,2)+Sij(1,3)+Sij(2,3)));
Kvrh = mean([Kvoigt,Kreuss]);
Gvoigt = (1/15)*(Cij(1,1)+Cij(2,2)+Cij(3,3)+3*(Cij(4,4)+Cij(5,5)+Cij(6,6))-Cij(1,2)-Cij(1,3)-Cij(2,3));
Greuss = 15/(4*(Sij(1,1)+Sij(2,2)+Sij(3,3)-Sij(1,2)-Sij(1,3)-Sij(2,3))+3*(Sij(4,4)+Sij(5,5)+Sij(6,6)));
Gvrh = mean([Gvoigt,Greuss]);

VpIso = ((Kvrh+(4/3)*Gvrh)/(rho))^0.5*10;
VsIso = (Gvrh/(rho))^0.5*10;
VpVsIso = VpIso/VsIso;

end
%% Plotting
%% Colormaps
cmap1 = [0 1 0.25 0;0.49 1 1 0;0.5 0.5 1 1;0.51 0 1 1;1 0.5 0 1];
cmapA = [0 1 1 1;0.0025 1 0.25 0;0.49 1 1 0;0.5 0.5 1 1;0.51 0 1 1;1 0.5 0 1];
cmap = zeros(64,3);
cmap2 = zeros(64,3);
CMPTS = linspace(0,1,64)';
CMPTS2 = linspace(0,1,64)';

for i = 1:3
    cmap(:,i) = interp1(cmap1(:,1),cmap1(:,i+1),CMPTS);
    cmap2(:,i) = interp1(cmapA(:,1),cmapA(:,i+1),CMPTS2);
end
%% Quiver plot for Vs1 polarization
% Num Qs defines how many degrees between quivers, if equal to Res it will
% put a quiver at every direction. If you want every other direction,
% multiply Res by 2, every third would be Res*3 etc.
NumQs = Res; 
VpPol1 = VpPolxyz(1:NumQs:end,1:NumQs:end,:);
Vs1Pol1 = Vs1Polxyz(1:NumQs:end,1:NumQs:end,:);
Xxpol = Xx(1:NumQs:end,1:NumQs:end);
Yypol = Yy(1:NumQs:end,1:NumQs:end);
Zzpol = Zz(1:NumQs:end,1:NumQs:end);
Vs1PolX = Vs1Pol1(:,:,1);
Vs1PolY = Vs1Pol1(:,:,2);
Vs1PolZ = Vs1Pol1(:,:,3);

% Patch shape to cover up the upper hemisphere part of the projection
P1 = [cos(-pi:pi/180:pi)',sin(-pi:pi/180:pi)';[-1.05,-1.05,1.05,1.05,-1.05]',[0,1.05,1.05,-1.05,-1.05]'];

% Turning x,y,z and velp, vels1, and vels2 matrices into column vectors
Uxyz = [reshape(Xx,[size(Xx,1)*size(Xx,2),1]),reshape(Yy,[size(Yy,1)*size(Yy,2),1]),reshape(Zz,[size(Zz,1)*size(Zz,2),1])];
Vp = reshape(velp,[size(velp,1)*size(velp,2),1]);
Vs1 = reshape(vels1,[size(vels1,1)*size(vels1,2),1]);
Vs2 = reshape(vels2,[size(vels2,1)*size(vels2,2),1]);

% Need to get rid of Z = 1 because makes an inf value for the projection
I = find(Uxyz(:,3)<1);
Uxyz = Uxyz(I,:);
Vp = Vp(I);
Vs1 = Vs1(I);
Vs2 = Vs2(I);

h = figure;
%set(h, 'Color', 'w', 'Name', FigName,'NumberTitle','off')
colormap(cmap2)
[Xtest,Ytest] = meshgrid(-1:0.01:1,-1:0.01:1);

Test = TriScatteredInterp(Uxyz(:,1).*real(sqrt(1./(1-Uxyz(:,3)))),Uxyz(:,2).*real(sqrt(1./(1-Uxyz(:,3)))),1./Vs2-1./Vs1);
subplot(2,3,3)
contourf(Xtest,Ytest,Test(Xtest,Ytest),10,'EdgeColor','none');
caxis([min(1./Vs2-1./Vs1) max(1./Vs2-1./Vs1)])
axis equal    
axis([-1.01 1.01 -1.01 1.01])
hold on

quiver(Xxpol.*real(sqrt(1./(1-Zzpol))),Yypol.*real(sqrt(1./(1-Zzpol))),...
     Vs1PolX,Vs1PolY,'k','ShowArrowHead','off')

patch(P1(:,1),P1(:,2),'w','EdgeColor','none')

plot(cos(-pi:pi/180:pi),sin(-pi:pi/180:pi),'k','LineWidth',1)
axis off
colorbar

Title = {[' Vs1 Polarization & splitting time (s/km)'];['Max: ',num2str(max(1./Vs2-1./Vs1)),' Min: ',num2str(min(1./Vs2-1./Vs1))];...
    ['AVsMax = ',num2str(max(Vs1-Vs2)/((max(Vs1)+min(Vs1))/2)*100)]};
title(Title)

% Vp plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
[Xtest,Ytest] = meshgrid(-1:0.01:1,-1:0.01:1);

Test = TriScatteredInterp(Uxyz(:,1).*real(sqrt(1./(1-Uxyz(:,3)))),Uxyz(:,2).*real(sqrt(1./(1-Uxyz(:,3)))),Vp);

subplot(2,3,1)
contourf(Xtest,Ytest,Test(Xtest,Ytest),10,'EdgeColor','none');
caxis([min(Vp) max(Vp)])
axis equal    
axis([-1.01 1.01 -1.01 1.01])

hold on
patch(P1(:,1),P1(:,2),'w','EdgeColor','none')
plot(cos(-pi:pi/180:pi),sin(-pi:pi/180:pi),'k','LineWidth',1)
axis off
colorbar

Title = {[' Vp (km/s): Max: ',num2str(max(Vp))];[' Min: ',num2str(min(Vp)),' AVp: ',num2str((max(Vp)-min(Vp))/((max(Vp)+min(Vp))/2)*100)]};
title(Title)

% Vs1 plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
[Xtest,Ytest] = meshgrid(-1:0.01:1,-1:0.01:1);

Test = TriScatteredInterp(Uxyz(:,1).*real(sqrt(1./(1-Uxyz(:,3)))),Uxyz(:,2).*real(sqrt(1./(1-Uxyz(:,3)))),Vs1);

subplot(2,3,4)
contourf(Xtest,Ytest,Test(Xtest,Ytest),10,'EdgeColor','none');
caxis([min(Vs1) max(Vs1)])
axis equal    
axis([-1.01 1.01 -1.01 1.01])

hold on
patch(P1(:,1),P1(:,2),'w','EdgeColor','none')
plot(cos(-pi:pi/180:pi),sin(-pi:pi/180:pi),'k','LineWidth',1)
axis off
colorbar

Title = {[' Vs1 (km/s): Max: ',num2str(max(Vs1))];[' Min: ',num2str(min(Vs1)),' AVs1: ',num2str((max(Vs1)-min(Vs1))/((max(Vs1)+min(Vs1))/2)*100)]};
title(Title)

% Vs2 plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
[Xtest,Ytest] = meshgrid(-1:0.01:1,-1:0.01:1);

Test = TriScatteredInterp(Uxyz(:,1).*real(sqrt(1./(1-Uxyz(:,3)))),Uxyz(:,2).*real(sqrt(1./(1-Uxyz(:,3)))),Vs2);

subplot(2,3,5)
contourf(Xtest,Ytest,Test(Xtest,Ytest),10,'EdgeColor','none');
caxis([min(Vs2) max(Vs2)])
axis equal    
axis([-1.01 1.01 -1.01 1.01])

hold on
patch(P1(:,1),P1(:,2),'w','EdgeColor','none')
plot(cos(-pi:pi/180:pi),sin(-pi:pi/180:pi),'k','LineWidth',1)
axis off
colorbar

Title = {[' Vs2 (km/s): Max: ',num2str(max(Vs2))];[' Min: ',num2str(min(Vs2)),' AVs2: ',num2str((max(Vs2)-min(Vs2))/((max(Vs2)+min(Vs2))/2)*100)]};
title(Title)

% Vp/Vs1 plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
[Xtest,Ytest] = meshgrid(-1:0.01:1,-1:0.01:1);

Test = TriScatteredInterp(Uxyz(:,1).*real(sqrt(1./(1-Uxyz(:,3)))),Uxyz(:,2).*real(sqrt(1./(1-Uxyz(:,3)))),Vp./Vs1);

subplot(2,3,2)
contourf(Xtest,Ytest,Test(Xtest,Ytest),10,'EdgeColor','none');
caxis([min(Vp./Vs1) max(Vp./Vs1)])
axis equal    
axis([-1.01 1.01 -1.01 1.01])

hold on
patch(P1(:,1),P1(:,2),'w','EdgeColor','none')
plot(cos(-pi:pi/180:pi),sin(-pi:pi/180:pi),'k','LineWidth',1)
axis off
colorbar

Title = {[' VpVs: Max: ',num2str(max(Vp./Vs1))];[' Min: ',num2str(min(Vp./Vs1))]};
title(Title)
