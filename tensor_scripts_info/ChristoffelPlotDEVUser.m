%% Modified Code from David Healy
% This code solves the Christoffel equation using direction cosines
% produced by the sphere function in Matlab. Default is every 3 degrees,
% (sphere(120)), but this can be changed for higher resolution by changing
% the Res variable here.
Res = 3;

%%
% Reads in data from xls file, with sample # in first column, then density in
% g/cm^3, then the 21 tensor components:
% Samp#,rho, C11, C12, C13, C14, C15, C16, C22, C23,...

%%% William Note 9/14/2020 - AllAggTensors.xslx is the input file needed to
%%% run this
%filename='user_input.xlsx';
filename=string(filename);
Data = xlsread(filename);
%% Choosing which tensors to plot
% 1 = Voigt, 2 = Reuss, 3 = VRH (only VRH for whole rock)
% Choose phase (0 = whole rock, 1 = quartz, 2 = plagioclase,...)
rockType = 1;

% change this to choose your sample
tensorN = 1;
FigName = ['Sample ',num2str(tensorN),' Phase = ',num2str(rockType)];
 
Tensors = Data((tensorN),4:24);

Rhos = Data((tensorN),3);
%Snum = Data(Min(VRH(Samps)),1);
%Samples = num2str(Snum);
%%

for s = 1:length(Rhos)
C = [Tensors(s,1) Tensors(s,2) Tensors(s,3) Tensors(s,4) Tensors(s,5) Tensors(s,6);...
    Tensors(s,2) Tensors(s,7) Tensors(s,8) Tensors(s,9) Tensors(s,10) Tensors(s,11);...
    Tensors(s,3) Tensors(s,8) Tensors(s,12) Tensors(s,13) Tensors(s,14) Tensors(s,15);...
    Tensors(s,4) Tensors(s,9) Tensors(s,13) Tensors(s,16) Tensors(s,17) Tensors(s,18);...
    Tensors(s,5) Tensors(s,10) Tensors(s,14) Tensors(s,17) Tensors(s,19) Tensors(s,20);...
    Tensors(s,6) Tensors(s,11) Tensors(s,15) Tensors(s,18) Tensors(s,20) Tensors(s,21)];

rho = Rhos(s);
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
set(h, 'Color', 'w', 'Name', FigName,'NumberTitle','off')
colormap(cmap2)
[Xtest,Ytest] = meshgrid(-1:0.01:1,-1:0.01:1);

Test = TriScatteredInterp(Uxyz(:,1).*real(sqrt(1./(1-Uxyz(:,3)))),Uxyz(:,2).*real(sqrt(1./(1-Uxyz(:,3)))),1./Vs2-1./Vs1);
%subplot(2,3,3)
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
%colorbar

Title = {[tensorN,' Vs1 Polarization & splitting time (s/km)'];['Max: ',num2str(max(1./Vs2-1./Vs1)),' Min: ',num2str(min(1./Vs2-1./Vs1))];...
    ['AVsMax = ',num2str(max(Vs1-Vs2)/((max(Vs1)+min(Vs1))/2)*100)]};
title(Title)

fig2plotly(gcf, 'offline', true, 'filename', 'quiver_plot')

%% Vp plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
[Xtest,Ytest] = meshgrid(-1:0.01:1,-1:0.01:1);

Test = TriScatteredInterp(Uxyz(:,1).*real(sqrt(1./(1-Uxyz(:,3)))),Uxyz(:,2).*real(sqrt(1./(1-Uxyz(:,3)))),Vp);

subplot(1,1,1)
contourf(Xtest,Ytest,Test(Xtest,Ytest),10,'EdgeColor','none');
caxis([min(Vp) max(Vp)])
axis equal    
axis([-1.01 1.01 -1.01 1.01])

hold on
patch(P1(:,1),P1(:,2),'w','EdgeColor','none')
plot(cos(-pi:pi/180:pi),sin(-pi:pi/180:pi),'k','LineWidth',1)
axis off
%colorbar

Title = {[tensorN,' Vp (km/s): Max: ',num2str(max(Vp))];[' Min: ',num2str(min(Vp)),' AVp: ',num2str((max(Vp)-min(Vp))/((max(Vp)+min(Vp))/2)*100)]};
title(Title)

fig2plotly(gcf, 'offline', true, 'filename', 'vp_plot')
%% Vs1 plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
[Xtest,Ytest] = meshgrid(-1:0.01:1,-1:0.01:1);

Test = TriScatteredInterp(Uxyz(:,1).*real(sqrt(1./(1-Uxyz(:,3)))),Uxyz(:,2).*real(sqrt(1./(1-Uxyz(:,3)))),Vs1);

%subplot(1,1,4)
contourf(Xtest,Ytest,Test(Xtest,Ytest),10,'EdgeColor','none');
caxis([min(Vs1) max(Vs1)])
axis equal    
axis([-1.01 1.01 -1.01 1.01])

hold on
patch(P1(:,1),P1(:,2),'w','EdgeColor','none')
plot(cos(-pi:pi/180:pi),sin(-pi:pi/180:pi),'k','LineWidth',1)
axis off
%colorbar

Title = {[tensorN,' Vs1 (km/s): Max: ',num2str(max(Vs1))];[' Min: ',num2str(min(Vs1)),' AVs1: ',num2str((max(Vs1)-min(Vs1))/((max(Vs1)+min(Vs1))/2)*100)]};
title(Title)

fig2plotly(gcf, 'offline', true, 'filename', 'vs1_plot')
%% Vs2 plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
[Xtest,Ytest] = meshgrid(-1:0.01:1,-1:0.01:1);

Test = TriScatteredInterp(Uxyz(:,1).*real(sqrt(1./(1-Uxyz(:,3)))),Uxyz(:,2).*real(sqrt(1./(1-Uxyz(:,3)))),Vs2);

%subplot(1,1,5)
contourf(Xtest,Ytest,Test(Xtest,Ytest),10,'EdgeColor','none');
caxis([min(Vs2) max(Vs2)])
axis equal    
axis([-1.01 1.01 -1.01 1.01])

hold on
patch(P1(:,1),P1(:,2),'w','EdgeColor','none')
plot(cos(-pi:pi/180:pi),sin(-pi:pi/180:pi),'k','LineWidth',1)
axis off
%colorbar

Title = {[tensorN,' Vs2 (km/s): Max: ',num2str(max(Vs2))];[' Min: ',num2str(min(Vs2)),' AVs2: ',num2str((max(Vs2)-min(Vs2))/((max(Vs2)+min(Vs2))/2)*100)]};
title(Title)

fig2plotly(gcf, 'offline', true, 'filename', 'vs2_plot')
%% Vp/Vs1 plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
[Xtest,Ytest] = meshgrid(-1:0.01:1,-1:0.01:1);

Test = TriScatteredInterp(Uxyz(:,1).*real(sqrt(1./(1-Uxyz(:,3)))),Uxyz(:,2).*real(sqrt(1./(1-Uxyz(:,3)))),Vp./Vs1);

%subplot(2,3,2)
contourf(Xtest,Ytest,Test(Xtest,Ytest),10,'EdgeColor','none');
caxis([min(Vp./Vs1) max(Vp./Vs1)])
axis equal    
axis([-1.01 1.01 -1.01 1.01])

hold on
patch(P1(:,1),P1(:,2),'w','EdgeColor','none')
plot(cos(-pi:pi/180:pi),sin(-pi:pi/180:pi),'k','LineWidth',1)
axis off
%colorbar

Title = {[tensorN,' VpVs: Max: ',num2str(max(Vp./Vs1))];[' Min: ',num2str(min(Vp./Vs1))]};
title(Title)
fig2plotly(gcf, 'offline', true, 'filename', 'vpvs1_plot')
%% Interactive plots
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

IntFigName = [FigName,' 3D Interactive'];
h = figure;
set(h, 'Color', 'w', 'Name', IntFigName,'NumberTitle','off')
colormap(cmap2)

%subplot(2,3,3)
surf(Xx,Yy,Zz,(1./vels2-1./vels1),'EdgeColor','flat')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
hold on
quiver3(Xxpol,Yypol,Zzpol,...
     Vs1PolX,Vs1PolY,Vs1PolZ,'k','ShowArrowHead','off')
% this patch is meant to be the foliation plane, right now it is a XZ
% plane, but we will probably want to have this orientation as an input
patch([1,1,-1,-1],[0,0,0,0],[1,-1,-1,1],[0.5 0.5 0.5])
% This dot is meant to be the lineation orientation, which is set to X
% right now
plot3([1,-1],[0,0],[0,0],'ko','MarkerFaceColor','k')
axis equal
axis off

%colorbar

Title = {[tensorN,' Vs1 Polarization & splitting time (s/km)'];['Max: ',num2str(max(1./Vs2-1./Vs1)),' Min: ',num2str(min(1./Vs2-1./Vs1))];...
    ['AVsMax = ',num2str(max(Vs1-Vs2)/((max(Vs1)+min(Vs1))/2)*100)]};
title(Title)

fig2plotly(gcf, 'offline', true, 'filename', 'vqvs1_plot')
%% Vp plot
%h = figure;
set(h, 'Color', 'w')
colormap(cmap)
%subplot(2,3,1)
surf(Xx,Yy,Zz,velp,'EdgeColor','flat')
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar
hold on
patch([1,1,-1,-1],[0,0,0,0],[1,-1,-1,1],[0.5 0.5 0.5])
plot3([1,-1],[0,0],[0,0],'ko','MarkerFaceColor','k')
axis equal
axis off

Title = {[tensorN,' Vp (km/s): Max: ',num2str(max(Vp))];[' Min: ',num2str(min(Vp)),' AVp: ',num2str((max(Vp)-min(Vp))/((max(Vp)+min(Vp))/2)*100)]};
title(Title)
%% Vs1 plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
%subplot(2,3,4)
surf(Xx,Yy,Zz,vels1,'EdgeColor','flat')
xlabel('X')
ylabel('Y')
zlabel('Z')
%colorbar
hold on
patch([1,1,-1,-1],[0,0,0,0],[1,-1,-1,1],[0.5 0.5 0.5])
plot3([1,-1],[0,0],[0,0],'ko','MarkerFaceColor','k')
axis equal
axis off

Title = {[tensorN,' Vs1 (km/s): Max: ',num2str(max(Vs1))];[' Min: ',num2str(min(Vs1)),' AVs1: ',num2str((max(Vs1)-min(Vs1))/((max(Vs1)+min(Vs1))/2)*100)]};
title(Title)
%% Vs2 plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
%subplot(2,3,5)
surf(Xx,Yy,Zz,vels2,'EdgeColor','flat')
xlabel('X')
ylabel('Y')
zlabel('Z')
%colorbar
hold on
patch([1,1,-1,-1],[0,0,0,0],[1,-1,-1,1],[0.5 0.5 0.5])
plot3([1,-1],[0,0],[0,0],'ko','MarkerFaceColor','k')
axis equal
axis off

Title = {[tensorN,' Vs2 (km/s): Max: ',num2str(max(Vs2))];[' Min: ',num2str(min(Vs2)),' AVs2: ',num2str((max(Vs2)-min(Vs2))/((max(Vs2)+min(Vs2))/2)*100)]};
title(Title)
%% Vp/Vs1 plot
%h = figure;
%set(h, 'Color', 'w')
colormap(cmap)
%subplot(2,3,2)
surf(Xx,Yy,Zz,velp./vels1,'EdgeColor','flat')
xlabel('X')
ylabel('Y')
zlabel('Z')
%colorbar
hold on
patch([1,1,-1,-1],[0,0,0,0],[1,-1,-1,1],[0.5 0.5 0.5])
plot3([1,-1],[0,0],[0,0],'ko','MarkerFaceColor','k')
axis equal
axis off

Title = {[tensorN,' VpVs: Max: ',num2str(max(Vp./Vs1))];[' Min: ',num2str(min(Vp./Vs1))]};
title(Title)
%% Backazimuthal plots for horizontal propagation

VpPolX = VpPolxyz(:,:,1);
VpPolY = VpPolxyz(:,:,2);
VpPolZ = VpPolxyz(:,:,3);
Vs1PolX = Vs1Polxyz(:,:,1);
Vs1PolY = Vs1Polxyz(:,:,2);
Vs1PolZ = Vs1Polxyz(:,:,3);
Vs2PolX = Vs2Polxyz(:,:,1);
Vs2PolY = Vs2Polxyz(:,:,2);
Vs2PolZ = Vs2Polxyz(:,:,3);

% Reshaping to column vectors
Uxyz = [reshape(Xx,size(Vs1PolZ,1)*size(Vs1PolZ,2),1),reshape(Yy,size(Vs1PolZ,1)*size(Vs1PolZ,2),1),reshape(Zz,size(Vs1PolZ,1)*size(Vs1PolZ,2),1)];

VS1baz = reshape(vels1,size(Vs1PolZ,1)*size(Vs1PolZ,2),1);
VS2baz = reshape(vels2,size(Vs1PolZ,1)*size(Vs1PolZ,2),1);
VPbaz = reshape(velp,size(Vs1PolZ,1)*size(Vs1PolZ,2),1);

VpPolZbaz = reshape(VpPolZ,size(Vs1PolZ,1)*size(Vs1PolZ,2),1);
Vs1PolZbaz = reshape(Vs1PolZ,size(Vs1PolZ,1)*size(Vs1PolZ,2),1);
Vs2PolZbaz = reshape(Vs2PolZ,size(Vs1PolZ,1)*size(Vs1PolZ,2),1);
VpPolXbaz = reshape(VpPolX,size(Vs1PolZ,1)*size(Vs1PolZ,2),1);
Vs1PolXbaz = reshape(Vs1PolX,size(Vs1PolX,1)*size(Vs1PolX,2),1);
Vs2PolXbaz = reshape(Vs2PolX,size(Vs1PolX,1)*size(Vs1PolX,2),1);
VpPolYbaz = reshape(VpPolY,size(Vs1PolZ,1)*size(Vs1PolZ,2),1);
Vs1PolYbaz = reshape(Vs1PolY,size(Vs1PolY,1)*size(Vs1PolY,2),1);
Vs2PolYbaz = reshape(Vs2PolY,size(Vs1PolY,1)*size(Vs1PolY,2),1);

% Take only horizontally propagating waves, or closest to if none at zero
HorzInt = 0.001;
Horzxyz = Uxyz(abs(Uxyz(:,3))<HorzInt,:);
Lxy = sqrt(Horzxyz(:,1).^2+Horzxyz(:,2).^2);
VS1baz = VS1baz(abs(Uxyz(:,3))<HorzInt);
VS2baz = VS2baz(abs(Uxyz(:,3))<HorzInt);
VPbaz = VPbaz(abs(Uxyz(:,3))<HorzInt);

Vs1PolZbaz = Vs1PolZbaz(abs(Uxyz(:,3))<HorzInt);
Vs2PolZbaz = Vs2PolZbaz(abs(Uxyz(:,3))<HorzInt);
Vs1PolXbaz = Vs1PolXbaz(abs(Uxyz(:,3))<HorzInt);
Vs2PolXbaz = Vs2PolXbaz(abs(Uxyz(:,3))<HorzInt);
Vs1PolYbaz = Vs1PolYbaz(abs(Uxyz(:,3))<HorzInt);
Vs2PolYbaz = Vs2PolYbaz(abs(Uxyz(:,3))<HorzInt);
VpPolXbaz = VpPolXbaz(abs(Uxyz(:,3))<HorzInt);
VpPolYbaz = VpPolYbaz(abs(Uxyz(:,3))<HorzInt);
VpPolZbaz = VpPolZbaz(abs(Uxyz(:,3))<HorzInt);

% calculate backazimuthal angles with y as North
BackAz = zeros(length(VPbaz),1);
for i = 1:length(BackAz)
    if Horzxyz(i,2)>=0
        BackAz(i) = asind(Horzxyz(i,1)/Lxy(i));
    elseif Horzxyz(i,2)<0
        BackAz(i) = 180 - asind(Horzxyz(i,1)/Lxy(i));
    end
end
BackAz(BackAz<0) = BackAz(BackAz<0) + 360;

% Red-blue colormap
RB = [0 1 0 0;0.5 1 1 1;1 0 0 1];
cmapRB = zeros(16,3);
CMPTSRB = linspace(0,1,16)';
for i = 1:3
    cmapRB(:,i) = interp1(RB(:,1),RB(:,i+1),CMPTSRB);
end

figure
set(gcf,'Color','w')
colormap(cmapRB)

% Inclination angles from horizontal (if this is 90 wave is vertically
% polarized
Vs1Inc = abs(atand(Vs1PolZbaz./(Vs1PolXbaz.^2+Vs1PolYbaz.^2).^0.5));
Vs2Inc = abs(atand(Vs2PolZbaz./(Vs2PolXbaz.^2+Vs2PolYbaz.^2).^0.5));

scatter(BackAz,VS1baz,10,abs(Vs1Inc),'filled')
hold on
scatter(BackAz,VS2baz,10,abs(Vs2Inc),'filled')
caxis([0 90])
xlabel('Backazimuth (?)')
ylabel('Vs (km/s)')
%colorbar
Title = {['Vs with Backazimuth and inclination from horizontal']};
title(Title)

fig2plotly(gcf, 'offline', true, 'filename', 'backazimuthal_plot')

%% Radial and Transvers components?
% Rotate polarization vectors by backazimuth around Z. I think this should
% make it like looking directly down the backazimuth...
R_PV1 = zeros(length(BackAz),3);
R_PV2 = zeros(length(BackAz),3);
R_PVP = zeros(length(BackAz),3);

for i = 1:length(BackAz)
    % Rotation matrix
    R = [cosd(BackAz(i)),-sind(BackAz(i)),0;...
        sind(BackAz(i)),cosd(BackAz(i)),0;...
        0,0,1];
    
    % Polarization vector
    PV1 = [Vs1PolXbaz(i),Vs1PolYbaz(i),Vs1PolZbaz(i)]';
    PV2 = [Vs2PolXbaz(i),Vs2PolYbaz(i),Vs2PolZbaz(i)]';
    PVP = [VpPolXbaz(i),VpPolYbaz(i),VpPolZbaz(i)]';
    
    R_PV1(i,:) = (R*PV1)';
    R_PV2(i,:) = (R*PV2)';
    R_PVP(i,:) = (R*PVP)';
end
%
figure
set(gcf,'Color','w')
VS1T = R_PV1(:,1);
VS1R = R_PV1(:,2);
VS1V = R_PV1(:,3);

VS2T = R_PV2(:,1);
VS2R = R_PV2(:,2);
VS2V = R_PV2(:,3);

VPT = R_PVP(:,1);
VPR = R_PVP(:,2);
VPV = R_PVP(:,3);

colormap(cmapRB)
%subplot(2,1,1),scatter(BackAz,VS1baz,10,abs(VS1V)-abs(VS1T),'filled')
axis([0,360,min(Vs2)-0.05,max(Vs1)+0.05])
caxis([-1 1])
hold on
scatter(BackAz,VS2baz,10,abs(VS2V)-abs(VS2T),'filled')
title('Vs (km/s) colored by polarization components (V-T)')
%colorbar
[I,BAZsort] = sort(BackAz);

%subplot(2,1,2),plot(BackAz(BAZsort),abs(VS1R(BAZsort)),'k','LineWidth',4)
% colorbar just so both plots are the same width
%colorbar
axis([0,360,0,1])
hold on
plot(BackAz(BAZsort),abs(VS1T(BAZsort)),'g','LineWidth',4)
plot(BackAz(BAZsort),abs(VS1V(BAZsort)),'b','LineWidth',4)
plot(BackAz(BAZsort),abs(VS2R(BAZsort)),'K:','LineWidth',2)
plot(BackAz(BAZsort),abs(VS2T(BAZsort)),'g:','LineWidth',2)
plot(BackAz(BAZsort),abs(VS2V(BAZsort)),'b:','LineWidth',2)
legend({'VS1 R','VS1 T','VS1 V','VS2 R','VS2 T','VS2 V'},'Location','northeast')
title('Polarization components (R, V, T) for Vs1 and Vs2')

fig2plotly(gcf, 'offline', true, 'filename', 'radial_plot')
exit
