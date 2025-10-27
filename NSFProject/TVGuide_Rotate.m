%%
close all
clear all
AllTensors = xlsread('PelonaDB_XYZorientation.xlsx');

% As it is now, this rotates all the tensors the same, and stores them in RotatedCij
RotatedCij = zeros(length(AllTensors(:,1)),23);
for t = 1:length(AllTensors(:,1))
    
C = AllTensors(t,3:23);
C = [C(1:6);...
    C(2), C(7:11);...
    C(3), C(8), C(12:15);...
    C(4), C(9), C(13), C(16:18);...
    C(5), C(10), C(14), C(17), C(19:20);...
    C(6), C(11), C(15), C(18), C(20), C(21)];
rho = AllTensors(:,2);

%% Enter rotation values here
dip = 30; % rotates around Y
strike = 200; % rotates around Z after applying dip to put strike at correct azimuth
% Using right hand rule, trend MUST be between strike and strike+180
trend = 270; % rotates by rake (calculated below) around Z prior to applying strike and dip to put lineation at correct azimuthal trend (Y is North)
% NOTE: because Y is North, to apply no rotation at all, you need to enter
% 90 for the trend, and 0 for both strike and dip.

% If dip = 90 (vertical foliation), trend must equal strike, so plunge can
% be used to specify within plane rotation of lineation. This is ONLY used
% if dip = 90
plunge = 45; 

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
%%

for i = 1:length(dip)
    ya(i) = dip(i); % Rotates around top, around Y if XY orientation, DIP
    xa(i) = 0; % Rotates around right, around X ''
    za(i) = strike(i); % Rotates around center, around Z '', STRIKE
    yrake(i) = rake(i)-90; % This also rotates around center, before the dip rotation, to put lineation at the right rake. The -90 is because Y is North.
end

%% We could put in options to apply multiple rotations

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
        
        end
        RotatedCij(t,1) = AllTensors(t,1);
        RotatedCij(t,2) = rho(t);
        RotatedCij(t,3:8) = CR(1,:);
        RotatedCij(t,9:13) = CR(2,2:6);
        RotatedCij(t,14:17) = CR(3,3:6);
        RotatedCij(t,18:20) = CR(4,4:6);
        RotatedCij(t,21:22) = CR(5,5:6);
        RotatedCij(t,23) = CR(6,6);
        

end
