
AllTensors = xlsread('PelonaDB_XYZorientation.xlsx');

% Select tensors to average
N = [69,70,83];
Cijs = zeros(length(N),21);
Rhos = zeros(length(N),1);
for i = 1:length(N)
    Cijs(i,:) = AllTensors(AllTensors(:,1)==N(i),3:23);
    Rhos(i) = AllTensors(AllTensors(:,1)==N(i),2);
end

% Specify weights
Equal = 0; % If 1, then weights will all be equal, otherwise must specify a weight for each tensor
% These are the weights if specifying
wt = [2, 1, 1];
wt = wt./sum(wt); % Makes sure weights add to 1

if Equal == 1
    wt = ones(length(N),1)*(1/length(N));
else
end


%% Averaging V is for Voigt (using Cij), R is for Reuss (using inverse of Cij)
TensorsV = cell(length(N),1);
TensorsR = cell(length(N),1);
V_ave = zeros(6,6);
R_ave = zeros(6,6);

for i = 1:length(N)
    C = Cijs(i,1:end);
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
    PV = zeros(length(N),1);
    PR = zeros(length(N),1);
    for j = 1:length(N)
        PV(j) = wt(j)*TensorsV{j}(i);
        PR(j) = wt(j)*TensorsR{j}(i);
    end
    V_ave(i) = sum(PV);
    R_ave(i) = sum(PR);
end
% VRH average is the average of V and R stiffness tensors
VRH = (V_ave + inv(R_ave))/2;