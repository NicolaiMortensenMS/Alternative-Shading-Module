function output = Klems_testing()
% This script compares full IGU BSDFs from WINDOW with IGUs built up from
% single layer BSDFs from WINDOW 7, combined using the Klems method.
%
% Authors
% Nicolai Mortensen, MicroShade A/S
%
% Last Updated:
% 26-05-2020

%% Defining Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reading BSDF

% BSDF matrices have been generated in LBNL WINDOW for glazing layers and
% IGUs. These have been saved in the folder BSDF -> glasses as .mat files.
BSDF_filePath = fullfile(getExecutableFolder(), 'BSDF','glasses'); % Folder containing BSDFs for glasses
fileDirectories = dir(BSDF_filePath);


% Saving all files in BSDF struct
k = 0;
for i = 1:length(fileDirectories)
    if fileDirectories(i).isdir == 0 % tests for non-directory files.
        k = k+1;
        fileNames{k} = fileDirectories(i).name; % saves file names.
        BSDF_fileName = fullfile(BSDF_filePath, char(fileNames(k)));
        BSDF(k) = struct2cell(load(BSDF_fileName)); % loads BSDFs with file names.
    end
end

% Subdividing BSDF into Klems matrices (145x145)
for i = 1:length(BSDF)
    file = BSDF{1,i};
    k = [];
    for j = 1:length(file)
        if isnan(file(j))
            k = [k,j];
            if length(k) > 1
                BSDF{length(k),i} = file(k(end-1)+1:k(end)-1,:);
            end
        end
    end
    BSDF{length(k),i} = fileNames{i};
end

% Defining parameters of Klems matrices (T,R)

% Solar Parameters
for i = 1:length(BSDF)
    % Front transmittance
    Tf_sol{1,i} = BSDF{end,i};
    Tf_sol{2,i} = BSDF{2,i};
    % Back transmittance
    Tb_sol{1,i} = BSDF{end,i};
    Tb_sol{2,i} = BSDF{3,i};
    % Front reflectance
    Rf_sol{1,i} = BSDF{end,i};
    Rf_sol{2,i} = BSDF{4,i};
    % Back reflectance
    Rb_sol{1,i} = BSDF{end,i};
    Rb_sol{2,i} = BSDF{5,i};
end

% Visible Parameters
for i = 1:length(BSDF)
    % Front transmittance
    Tf_vis{1,i} = BSDF{end,i};
    Tf_vis{2,i} = BSDF{6,i};
    % Back transmittance
    Tb_vis{1,i} = BSDF{end,i};
    Tb_vis{2,i} = BSDF{7,i};
    % Front reflectance
    Rf_vis{1,i} = BSDF{end,i};
    Rf_vis{2,i} = BSDF{8,i};
    % Back reflectance
    Rb_vis{1,i} = BSDF{end,i};
    Rb_vis{2,i} = BSDF{9,i};
end

% Locating glass files
% Clear
PC4_loc = find(ismember(fileNames, 'PC4.mat'));
PC4_14_PC4_loc = find(ismember(fileNames, 'PC4-14-PC4.mat'));
PC4_14_PC4_14_PC4_loc = find(ismember(fileNames, 'PC4-14-PC4-14-PC4.mat'));
% PXN
PXN_PC4_loc = find(ismember(fileNames, '(PXN)PC4.mat'));
PC4_14_PXN_PC4_loc = find(ismember(fileNames, 'PC4-14-(PXN)PC4.mat'));
PC4_14_PXN_PC4_14_PXN_PC4_loc = find(ismember(fileNames, 'PC4-14-(PXN)PC4-14-(PXN)PC4.mat'));
% PO
PO_PC4_loc = find(ismember(fileNames, '(PO)PC4.mat'));
PC4_14_PO_PC4_loc = find(ismember(fileNames, 'PC4-14-(PO)PC4.mat'));
PC4_14_PO_PC4_14_PO_PC4_loc = find(ismember(fileNames, 'PC4-14-(PO)PC4-14-(PO)PC4.mat'));
% Solar Control
PC6_CLX6028_loc = find(ismember(fileNames, 'PC6(CLX6028).mat'));
PC6_CLX6028_14_PC4_loc = find(ismember(fileNames, 'PC6(CLX6028)-14-PC4.mat'));
PC6_CLX6028_14_PC4_14_PXN_PC4_loc = find(ismember(fileNames, 'PC6(CLX6028)-14-PC4-14-(PXN)PC4.mat'));
PC6_CLX6028_14_PC4_14_PO_PC4_loc = find(ismember(fileNames, 'PC6(CLX6028)-14-PC4-14-(PO)PC4.mat'));

%% Identity matrix (I)
I = eye(145);

%% Diagonal propagation matrix (D)

% The diagonal propagation matrix is defined in [Bueno2017], and is based
% on the minimum, maximum and mean polar angles of the 145 patches in the
% Klems BSDF structure. 

polar_min = [0,5,15,25,35,45,55,65,75];
polar_max = [polar_min(2:end),90];
polar_mean = (polar_min+polar_max)/2;
azimuth_divisions = [1,8,16,20,24,24,24,16,12];

N = [];
theta_min = [];
theta_max = [];
theta_mean = [];
for i = 1: length(azimuth_divisions)
    azimuth_divisions_cumulative(i) = sum(azimuth_divisions(1:i));
    N = [N,repelem(azimuth_divisions(i),azimuth_divisions(i))];
    theta_min = [theta_min,repelem(polar_min(i),azimuth_divisions(i))];
    theta_max = [theta_max,repelem(polar_max(i),azimuth_divisions(i))];
    theta_mean = [theta_mean,repelem(polar_mean(i),azimuth_divisions(i))];
end

D = zeros(145);
for i = 1:length(D)
    D(i,i) = 2*pi/N(i)*(cosd(theta_min(i))-cosd(theta_max(i)))*cosd(theta_mean(i));
end

%% Testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Glass

% Glass layers
% Layer 1 - PC4
T1_f = Tf_sol{end,PC4_loc};
T1_b = Tb_sol{end,PC4_loc};
R1_f = Rf_sol{end,PC4_loc};
R1_b = Rb_sol{end,PC4_loc};
% Layer 2 - PC4
T2_f = T1_f;
T2_b = T1_b;
R2_f = R1_f;
R2_b = R1_b;
% Layer 3 - PC4:
T3_f = T1_f;
T3_b = T1_b;
R3_f = R1_f;
R3_b = R1_b;

% 2-Layer Klems combination
T12_f = T2_f*inv(I-D*R1_b*D*R2_f)*D*T1_f;
R12_f = R1_f+T1_b*inv(I-D*R2_f*D*R1_b)*D*R2_f*D*T1_f;
T12_b = T1_b*inv(I-D*R2_f*D*R1_b)*D*T2_b;
R12_b = R2_b+T2_f*inv(I-D*R1_b*D*R2_f)*D*R1_b*D*T2_b;

% 2-layer pre-defined glazing system (ref)
T12_f_ref = Tf_sol{end,PC4_14_PC4_loc};
R12_f_ref = Rf_sol{end,PC4_14_PC4_loc};
T12_b_ref = Tb_sol{end,PC4_14_PC4_loc};
R12_b_ref = Rb_sol{end,PC4_14_PC4_loc};

% 3-Layer Klems combination
T123_f = T3_f*inv(I-D*R12_b*D*R3_f)*D*T12_f;
R123_f = R12_f+T12_b*inv(I-D*R3_f*D*R12_b)*D*R3_f*D*T12_f;
T123_b = T12_b*inv(I-D*R3_f*D*R12_b)*D*T3_b;
R123_b = R3_b+T3_f*inv(I-D*R12_b*D*R3_f)*D*R12_b*D*T3_b;

% 3-layer pre-defined glazing system (ref)
T123_f_ref = Tf_sol{end,PC4_14_PC4_14_PC4_loc};
R123_f_ref = Rf_sol{end,PC4_14_PC4_14_PC4_loc};
T123_b_ref = Tb_sol{end,PC4_14_PC4_14_PC4_loc};
R123_b_ref = Rb_sol{end,PC4_14_PC4_14_PC4_loc};

% Calculating relative differences
for i = 1:145
    dT12_f(i) = 100*abs(T12_f(i,i)-T12_f_ref(i,i))/((T12_f(i,i)+T12_f_ref(i,i))/2);
    dT12_b(i) = 100*abs(T12_b(i,i)-T12_b_ref(i,i))/((T12_b(i,i)+T12_b_ref(i,i))/2);
    dR12_f(i) = 100*abs(R12_f(i,i)-R12_f_ref(i,i))/((R12_f(i,i)+R12_f_ref(i,i))/2);
    dR12_b(i) = 100*abs(R12_b(i,i)-R12_b_ref(i,i))/((R12_b(i,i)+R12_b_ref(i,i))/2);
    dT123_f(i) = 100*abs(T123_f(i,i)-T123_f_ref(i,i))/((T123_f(i,i)+T123_f_ref(i,i))/2);
    dT123_b(i) = 100*abs(T123_b(i,i)-T123_b_ref(i,i))/((T123_b(i,i)+T123_b_ref(i,i))/2);
    dR123_f(i) = 100*abs(R123_f(i,i)-R123_f_ref(i,i))/((R123_f(i,i)+R123_f_ref(i,i))/2);
    dR123_b(i) = 100*abs(R123_b(i,i)-R123_b_ref(i,i))/((R123_b(i,i)+R123_b_ref(i,i))/2);
end

% Plotting relative differences

% Front transmitance - 2-layer
subplot(4,2,1)
x_axis = [1:1:145];
y_axis = dT12_f;
plot(x_axis,y_axis)
title('Front Transmittance - PC4-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back transmitance - 2-layer
subplot(4,2,3)
x_axis = [1:1:145];
y_axis = dT12_b;
plot(x_axis,y_axis)
title('Back Transmittance - PC4-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front reflectance - 2-layer
subplot(4,2,5)
x_axis = [1:1:145];
y_axis = dR12_f;
plot(x_axis,y_axis)
title('Front Reflectance - PC4-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back reflectance - 2-layer
subplot(4,2,7)
x_axis = [1:1:145];
y_axis = dR12_b;
plot(x_axis,y_axis)
title('Back Reflectance - PC4-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front transmitance - 3-layer
subplot(4,2,2)
x_axis = [1:1:145];
y_axis = dT123_f;
plot(x_axis,y_axis)
title('Front Transmittance - PC4-14-PC4-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back transmitance - 3-layer
subplot(4,2,4)
x_axis = [1:1:145];
y_axis = dT123_b;
plot(x_axis,y_axis)
title('Back Transmittance - PC4-14-PC4-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front reflectance - 3-layer
subplot(4,2,6)
x_axis = [1:1:145];
y_axis = dR123_f;
plot(x_axis,y_axis)
title('Front Reflectance - PC4-14-PC4-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back reflectance - 3-layer
subplot(4,2,8)
x_axis = [1:1:145];
y_axis = dR123_b;
plot(x_axis,y_axis)
title('Back Reflectance - PC4-14-PC4-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

%% PXN:

% Glass layers
% Layer 1 - PC4
T1_f = Tf_sol{end,PC4_loc};
T1_b = Tb_sol{end,PC4_loc};
R1_f = Rf_sol{end,PC4_loc};
R1_b = Rb_sol{end,PC4_loc};
% Layer 2 - (PXN)PC4
T2_f = Tf_sol{end,PXN_PC4_loc};
T2_b = Tb_sol{end,PXN_PC4_loc};
R2_f = Rf_sol{end,PXN_PC4_loc};
R2_b = Rb_sol{end,PXN_PC4_loc};
% Layer 3 - (PXN)PC4:
T3_f = T2_f;
T3_b = T2_b;
R3_f = R2_f;
R3_b = R2_b;

% 2-Layer Klems combination
T12_f = T2_f*inv(I-D*R1_b*D*R2_f)*D*T1_f;
R12_f = R1_f+T1_b*inv(I-D*R2_f*D*R1_b)*D*R2_f*D*T1_f;
T12_b = T1_b*inv(I-D*R2_f*D*R1_b)*D*T2_b;
R12_b = R2_b+T2_f*inv(I-D*R1_b*D*R2_f)*D*R1_b*D*T2_b;

% 2-layer pre-defined glazing system (ref)
T12_f_ref = Tf_sol{end,PC4_14_PXN_PC4_loc};
R12_f_ref = Rf_sol{end,PC4_14_PXN_PC4_loc};
T12_b_ref = Tb_sol{end,PC4_14_PXN_PC4_loc};
R12_b_ref = Rb_sol{end,PC4_14_PXN_PC4_loc};

% 3-Layer Klems combination
T123_f = T3_f*inv(I-D*R12_b*D*R3_f)*D*T12_f;
R123_f = R12_f+T12_b*inv(I-D*R3_f*D*R12_b)*D*R3_f*D*T12_f;
T123_b = T12_b*inv(I-D*R3_f*D*R12_b)*D*T3_b;
R123_b = R3_b+T3_f*inv(I-D*R12_b*D*R3_f)*D*R12_b*D*T3_b;

% 3-layer pre-defined glazing system (ref)
T123_f_ref = Tf_sol{end,PC4_14_PXN_PC4_14_PXN_PC4_loc};
R123_f_ref = Rf_sol{end,PC4_14_PXN_PC4_14_PXN_PC4_loc};
T123_b_ref = Tb_sol{end,PC4_14_PXN_PC4_14_PXN_PC4_loc};
R123_b_ref = Rb_sol{end,PC4_14_PXN_PC4_14_PXN_PC4_loc};

% Calculating relative differences
for i = 1:145
    dT12_f(i) = 100*abs(T12_f(i,i)-T12_f_ref(i,i))/((T12_f(i,i)+T12_f_ref(i,i))/2);
    dT12_b(i) = 100*abs(T12_b(i,i)-T12_b_ref(i,i))/((T12_b(i,i)+T12_b_ref(i,i))/2);
    dR12_f(i) = 100*abs(R12_f(i,i)-R12_f_ref(i,i))/((R12_f(i,i)+R12_f_ref(i,i))/2);
    dR12_b(i) = 100*abs(R12_b(i,i)-R12_b_ref(i,i))/((R12_b(i,i)+R12_b_ref(i,i))/2);
    dT123_f(i) = 100*abs(T123_f(i,i)-T123_f_ref(i,i))/((T123_f(i,i)+T123_f_ref(i,i))/2);
    dT123_b(i) = 100*abs(T123_b(i,i)-T123_b_ref(i,i))/((T123_b(i,i)+T123_b_ref(i,i))/2);
    dR123_f(i) = 100*abs(R123_f(i,i)-R123_f_ref(i,i))/((R123_f(i,i)+R123_f_ref(i,i))/2);
    dR123_b(i) = 100*abs(R123_b(i,i)-R123_b_ref(i,i))/((R123_b(i,i)+R123_b_ref(i,i))/2);
end

% Plotting relative differences

% Front transmitance - 2-layer
subplot(4,2,1)
x_axis = [1:1:145];
y_axis = dT12_f;
plot(x_axis,y_axis)
title('Front Transmittance - PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back transmitance - 2-layer
subplot(4,2,3)
x_axis = [1:1:145];
y_axis = dT12_b;
plot(x_axis,y_axis)
title('Back Transmittance - PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front reflectance - 2-layer
subplot(4,2,5)
x_axis = [1:1:145];
y_axis = dR12_f;
plot(x_axis,y_axis)
title('Front Reflectance - PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back reflectance - 2-layer
subplot(4,2,7)
x_axis = [1:1:145];
y_axis = dR12_b;
plot(x_axis,y_axis)
title('Back Reflectance - PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front transmitance - 3-layer
subplot(4,2,2)
x_axis = [1:1:145];
y_axis = dT123_f;
plot(x_axis,y_axis)
title('Front Transmittance - PC4-14-(PXN)PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back transmitance - 3-layer
subplot(4,2,4)
x_axis = [1:1:145];
y_axis = dT123_b;
plot(x_axis,y_axis)
title('Back Transmittance - PC4-14-(PXN)PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front reflectance - 3-layer
subplot(4,2,6)
x_axis = [1:1:145];
y_axis = dR123_f;
plot(x_axis,y_axis)
title('Front Reflectance - PC4-14-(PXN)PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back reflectance - 3-layer
subplot(4,2,8)
x_axis = [1:1:145];
y_axis = dR123_b;
plot(x_axis,y_axis)
title('Back Reflectance - PC4-14-(PXN)PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

%% PO:

% Glass layers
% Layer 1 - PC4
T1_f = Tf_sol{end,PC4_loc};
T1_b = Tb_sol{end,PC4_loc};
R1_f = Rf_sol{end,PC4_loc};
R1_b = Rb_sol{end,PC4_loc};
% Layer 2 - (PXN)PC4
T2_f = Tf_sol{end,PO_PC4_loc};
T2_b = Tb_sol{end,PO_PC4_loc};
R2_f = Rf_sol{end,PO_PC4_loc};
R2_b = Rb_sol{end,PO_PC4_loc};
% Layer 3 - (PXN)PC4:
T3_f = T2_f;
T3_b = T2_b;
R3_f = R2_f;
R3_b = R2_b;

% 2-Layer Klems combination
T12_f = T2_f*inv(I-D*R1_b*D*R2_f)*D*T1_f;
R12_f = R1_f+T1_b*inv(I-D*R2_f*D*R1_b)*D*R2_f*D*T1_f;
T12_b = T1_b*inv(I-D*R2_f*D*R1_b)*D*T2_b;
R12_b = R2_b+T2_f*inv(I-D*R1_b*D*R2_f)*D*R1_b*D*T2_b;

% 2-layer pre-defined glazing system (ref)
T12_f_ref = Tf_sol{end,PC4_14_PO_PC4_loc};
R12_f_ref = Rf_sol{end,PC4_14_PO_PC4_loc};
T12_b_ref = Tb_sol{end,PC4_14_PO_PC4_loc};
R12_b_ref = Rb_sol{end,PC4_14_PO_PC4_loc};

% 3-Layer Klems combination
T123_f = T3_f*inv(I-D*R12_b*D*R3_f)*D*T12_f;
R123_f = R12_f+T12_b*inv(I-D*R3_f*D*R12_b)*D*R3_f*D*T12_f;
T123_b = T12_b*inv(I-D*R3_f*D*R12_b)*D*T3_b;
R123_b = R3_b+T3_f*inv(I-D*R12_b*D*R3_f)*D*R12_b*D*T3_b;

% 3-layer pre-defined glazing system (ref)
T123_f_ref = Tf_sol{end,PC4_14_PO_PC4_14_PO_PC4_loc};
R123_f_ref = Rf_sol{end,PC4_14_PO_PC4_14_PO_PC4_loc};
T123_b_ref = Tb_sol{end,PC4_14_PO_PC4_14_PO_PC4_loc};
R123_b_ref = Rb_sol{end,PC4_14_PO_PC4_14_PO_PC4_loc};

% Calculating relative differences
for i = 1:145
    dT12_f(i) = 100*abs(T12_f(i,i)-T12_f_ref(i,i))/((T12_f(i,i)+T12_f_ref(i,i))/2);
    dT12_b(i) = 100*abs(T12_b(i,i)-T12_b_ref(i,i))/((T12_b(i,i)+T12_b_ref(i,i))/2);
    dR12_f(i) = 100*abs(R12_f(i,i)-R12_f_ref(i,i))/((R12_f(i,i)+R12_f_ref(i,i))/2);
    dR12_b(i) = 100*abs(R12_b(i,i)-R12_b_ref(i,i))/((R12_b(i,i)+R12_b_ref(i,i))/2);
    dT123_f(i) = 100*abs(T123_f(i,i)-T123_f_ref(i,i))/((T123_f(i,i)+T123_f_ref(i,i))/2);
    dT123_b(i) = 100*abs(T123_b(i,i)-T123_b_ref(i,i))/((T123_b(i,i)+T123_b_ref(i,i))/2);
    dR123_f(i) = 100*abs(R123_f(i,i)-R123_f_ref(i,i))/((R123_f(i,i)+R123_f_ref(i,i))/2);
    dR123_b(i) = 100*abs(R123_b(i,i)-R123_b_ref(i,i))/((R123_b(i,i)+R123_b_ref(i,i))/2);
end

% Plotting relative differences

% Front transmitance - 2-layer
subplot(4,2,1)
x_axis = [1:1:145];
y_axis = dT12_f;
plot(x_axis,y_axis)
title('Front Transmittance - PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back transmitance - 2-layer
subplot(4,2,3)
x_axis = [1:1:145];
y_axis = dT12_b;
plot(x_axis,y_axis)
title('Back Transmittance - PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front reflectance - 2-layer
subplot(4,2,5)
x_axis = [1:1:145];
y_axis = dR12_f;
plot(x_axis,y_axis)
title('Front Reflectance - PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back reflectance - 2-layer
subplot(4,2,7)
x_axis = [1:1:145];
y_axis = dR12_b;
plot(x_axis,y_axis)
title('Back Reflectance - PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front transmitance - 3-layer
subplot(4,2,2)
x_axis = [1:1:145];
y_axis = dT123_f;
plot(x_axis,y_axis)
title('Front Transmittance - PC4-14-(PO)PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back transmitance - 3-layer
subplot(4,2,4)
x_axis = [1:1:145];
y_axis = dT123_b;
plot(x_axis,y_axis)
title('Back Transmittance - PC4-14-(PO)PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front reflectance - 3-layer
subplot(4,2,6)
x_axis = [1:1:145];
y_axis = dR123_f;
plot(x_axis,y_axis)
title('Front Reflectance - PC4-14-(PO)PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back reflectance - 3-layer
subplot(4,2,8)
x_axis = [1:1:145];
y_axis = dR123_b;
plot(x_axis,y_axis)
title('Back Reflectance - PC4-14-(PO)PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

%% Solar Control + PXN:

% Glass layers
% Layer 1 - PC4
T1_f = Tf_sol{end,PC6_CLX6028_loc};
T1_b = Tb_sol{end,PC6_CLX6028_loc};
R1_f = Rf_sol{end,PC6_CLX6028_loc};
R1_b = Rb_sol{end,PC6_CLX6028_loc};
% Layer 2 - (PXN)PC4
T2_f = Tf_sol{end,PC4_loc};
T2_b = Tb_sol{end,PC4_loc};
R2_f = Rf_sol{end,PC4_loc};
R2_b = Rb_sol{end,PC4_loc};
% Layer 3 - (PXN)PC4:
T3_f = Tf_sol{end,PXN_PC4_loc};
T3_b = Tb_sol{end,PXN_PC4_loc};
R3_f = Rf_sol{end,PXN_PC4_loc};
R3_b = Rb_sol{end,PXN_PC4_loc};

% 2-Layer Klems combination
T12_f = T2_f*inv(I-D*R1_b*D*R2_f)*D*T1_f;
R12_f = R1_f+T1_b*inv(I-D*R2_f*D*R1_b)*D*R2_f*D*T1_f;
T12_b = T1_b*inv(I-D*R2_f*D*R1_b)*D*T2_b;
R12_b = R2_b+T2_f*inv(I-D*R1_b*D*R2_f)*D*R1_b*D*T2_b;

% 2-layer pre-defined glazing system (ref)
T12_f_ref = Tf_sol{end,PC6_CLX6028_14_PC4_loc};
R12_f_ref = Rf_sol{end,PC6_CLX6028_14_PC4_loc};
T12_b_ref = Tb_sol{end,PC6_CLX6028_14_PC4_loc};
R12_b_ref = Rb_sol{end,PC6_CLX6028_14_PC4_loc};

% 3-Layer Klems combination
T123_f = T3_f*inv(I-D*R12_b*D*R3_f)*D*T12_f;
R123_f = R12_f+T12_b*inv(I-D*R3_f*D*R12_b)*D*R3_f*D*T12_f;
T123_b = T12_b*inv(I-D*R3_f*D*R12_b)*D*T3_b;
R123_b = R3_b+T3_f*inv(I-D*R12_b*D*R3_f)*D*R12_b*D*T3_b;

% 3-layer pre-defined glazing system (ref)
T123_f_ref = Tf_sol{end,PC6_CLX6028_14_PC4_14_PXN_PC4_loc};
R123_f_ref = Rf_sol{end,PC6_CLX6028_14_PC4_14_PXN_PC4_loc};
T123_b_ref = Tb_sol{end,PC6_CLX6028_14_PC4_14_PXN_PC4_loc};
R123_b_ref = Rb_sol{end,PC6_CLX6028_14_PC4_14_PXN_PC4_loc};

% Calculating relative differences
for i = 1:145
    dT12_f(i) = 100*abs(T12_f(i,i)-T12_f_ref(i,i))/((T12_f(i,i)+T12_f_ref(i,i))/2);
    dT12_b(i) = 100*abs(T12_b(i,i)-T12_b_ref(i,i))/((T12_b(i,i)+T12_b_ref(i,i))/2);
    dR12_f(i) = 100*abs(R12_f(i,i)-R12_f_ref(i,i))/((R12_f(i,i)+R12_f_ref(i,i))/2);
    dR12_b(i) = 100*abs(R12_b(i,i)-R12_b_ref(i,i))/((R12_b(i,i)+R12_b_ref(i,i))/2);
    dT123_f(i) = 100*abs(T123_f(i,i)-T123_f_ref(i,i))/((T123_f(i,i)+T123_f_ref(i,i))/2);
    dT123_b(i) = 100*abs(T123_b(i,i)-T123_b_ref(i,i))/((T123_b(i,i)+T123_b_ref(i,i))/2);
    dR123_f(i) = 100*abs(R123_f(i,i)-R123_f_ref(i,i))/((R123_f(i,i)+R123_f_ref(i,i))/2);
    dR123_b(i) = 100*abs(R123_b(i,i)-R123_b_ref(i,i))/((R123_b(i,i)+R123_b_ref(i,i))/2);
end

% Plotting relative differences

% Front transmitance - 2-layer
subplot(4,2,1)
x_axis = [1:1:145];
y_axis = dT12_f;
plot(x_axis,y_axis)
title('Front Transmittance - PC6(CLX6028)-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back transmitance - 2-layer
subplot(4,2,3)
x_axis = [1:1:145];
y_axis = dT12_b;
plot(x_axis,y_axis)
title('Back Transmittance - PC6(CLX6028)-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front reflectance - 2-layer
subplot(4,2,5)
x_axis = [1:1:145];
y_axis = dR12_f;
plot(x_axis,y_axis)
title('Front Reflectance - PC6(CLX6028)-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back reflectance - 2-layer
subplot(4,2,7)
x_axis = [1:1:145];
y_axis = dR12_b;
plot(x_axis,y_axis)
title('Back Reflectance - PC6(CLX6028)-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front transmitance - 3-layer
subplot(4,2,2)
x_axis = [1:1:145];
y_axis = dT123_f;
plot(x_axis,y_axis)
title('Front Transmittance - PC6(CLX6028)-14-PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back transmitance - 3-layer
subplot(4,2,4)
x_axis = [1:1:145];
y_axis = dT123_b;
plot(x_axis,y_axis)
title('Back Transmittance - PC6(CLX6028)-14-PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front reflectance - 3-layer
subplot(4,2,6)
x_axis = [1:1:145];
y_axis = dR123_f;
plot(x_axis,y_axis)
title('Front Reflectance - PC6(CLX6028)-14-PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back reflectance - 3-layer
subplot(4,2,8)
x_axis = [1:1:145];
y_axis = dR123_b;
plot(x_axis,y_axis)
title('Back Reflectance - PC6(CLX6028)-14-PC4-14-(PXN)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

%% Solar Control + PO:

% Glass layers
% Layer 1 - PC4
T1_f = Tf_sol{end,PC6_CLX6028_loc};
T1_b = Tb_sol{end,PC6_CLX6028_loc};
R1_f = Rf_sol{end,PC6_CLX6028_loc};
R1_b = Rb_sol{end,PC6_CLX6028_loc};
% Layer 2 - (PXN)PC4
T2_f = Tf_sol{end,PC4_loc};
T2_b = Tb_sol{end,PC4_loc};
R2_f = Rf_sol{end,PC4_loc};
R2_b = Rb_sol{end,PC4_loc};
% Layer 3 - (PXN)PC4:
T3_f = Tf_sol{end,PO_PC4_loc};
T3_b = Tb_sol{end,PO_PC4_loc};
R3_f = Rf_sol{end,PO_PC4_loc};
R3_b = Rb_sol{end,PO_PC4_loc};

% 2-Layer Klems combination
T12_f = T2_f*inv(I-D*R1_b*D*R2_f)*D*T1_f;
R12_f = R1_f+T1_b*inv(I-D*R2_f*D*R1_b)*D*R2_f*D*T1_f;
T12_b = T1_b*inv(I-D*R2_f*D*R1_b)*D*T2_b;
R12_b = R2_b+T2_f*inv(I-D*R1_b*D*R2_f)*D*R1_b*D*T2_b;

% 2-layer pre-defined glazing system (ref)
T12_f_ref = Tf_sol{end,PC6_CLX6028_14_PC4_loc};
R12_f_ref = Rf_sol{end,PC6_CLX6028_14_PC4_loc};
T12_b_ref = Tb_sol{end,PC6_CLX6028_14_PC4_loc};
R12_b_ref = Rb_sol{end,PC6_CLX6028_14_PC4_loc};

% 3-Layer Klems combination
T123_f = T3_f*inv(I-D*R12_b*D*R3_f)*D*T12_f;
R123_f = R12_f+T12_b*inv(I-D*R3_f*D*R12_b)*D*R3_f*D*T12_f;
T123_b = T12_b*inv(I-D*R3_f*D*R12_b)*D*T3_b;
R123_b = R3_b+T3_f*inv(I-D*R12_b*D*R3_f)*D*R12_b*D*T3_b;

% 3-layer pre-defined glazing system (ref)
T123_f_ref = Tf_sol{end,PC6_CLX6028_14_PC4_14_PO_PC4_loc};
R123_f_ref = Rf_sol{end,PC6_CLX6028_14_PC4_14_PO_PC4_loc};
T123_b_ref = Tb_sol{end,PC6_CLX6028_14_PC4_14_PO_PC4_loc};
R123_b_ref = Rb_sol{end,PC6_CLX6028_14_PC4_14_PO_PC4_loc};

% Calculating relative differences
for i = 1:145
    dT12_f(i) = 100*abs(T12_f(i,i)-T12_f_ref(i,i))/((T12_f(i,i)+T12_f_ref(i,i))/2);
    dT12_b(i) = 100*abs(T12_b(i,i)-T12_b_ref(i,i))/((T12_b(i,i)+T12_b_ref(i,i))/2);
    dR12_f(i) = 100*abs(R12_f(i,i)-R12_f_ref(i,i))/((R12_f(i,i)+R12_f_ref(i,i))/2);
    dR12_b(i) = 100*abs(R12_b(i,i)-R12_b_ref(i,i))/((R12_b(i,i)+R12_b_ref(i,i))/2);
    dT123_f(i) = 100*abs(T123_f(i,i)-T123_f_ref(i,i))/((T123_f(i,i)+T123_f_ref(i,i))/2);
    dT123_b(i) = 100*abs(T123_b(i,i)-T123_b_ref(i,i))/((T123_b(i,i)+T123_b_ref(i,i))/2);
    dR123_f(i) = 100*abs(R123_f(i,i)-R123_f_ref(i,i))/((R123_f(i,i)+R123_f_ref(i,i))/2);
    dR123_b(i) = 100*abs(R123_b(i,i)-R123_b_ref(i,i))/((R123_b(i,i)+R123_b_ref(i,i))/2);
end

% Plotting relative differences

% Front transmitance - 2-layer
subplot(4,2,1)
x_axis = [1:1:145];
y_axis = dT12_f;
plot(x_axis,y_axis)
title('Front Transmittance - PC6(CLX6028)-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back transmitance - 2-layer
subplot(4,2,3)
x_axis = [1:1:145];
y_axis = dT12_b;
plot(x_axis,y_axis)
title('Back Transmittance - PC6(CLX6028)-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front reflectance - 2-layer
subplot(4,2,5)
x_axis = [1:1:145];
y_axis = dR12_f;
plot(x_axis,y_axis)
title('Front Reflectance - PC6(CLX6028)-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back reflectance - 2-layer
subplot(4,2,7)
x_axis = [1:1:145];
y_axis = dR12_b;
plot(x_axis,y_axis)
title('Back Reflectance - PC6(CLX6028)-14-PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front transmitance - 3-layer
subplot(4,2,2)
x_axis = [1:1:145];
y_axis = dT123_f;
plot(x_axis,y_axis)
title('Front Transmittance - PC6(CLX6028)-14-PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back transmitance - 3-layer
subplot(4,2,4)
x_axis = [1:1:145];
y_axis = dT123_b;
plot(x_axis,y_axis)
title('Back Transmittance - PC6(CLX6028)-14-PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Front reflectance - 3-layer
subplot(4,2,6)
x_axis = [1:1:145];
y_axis = dR123_f;
plot(x_axis,y_axis)
title('Front Reflectance - PC6(CLX6028)-14-PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

% Back reflectance - 3-layer
subplot(4,2,8)
x_axis = [1:1:145];
y_axis = dR123_b;
plot(x_axis,y_axis)
title('Back Reflectance - PC6(CLX6028)-14-PC4-14-(PO)PC4 - Relative Differences')
xlabel('BSDF Patch')
ylabel('Relative difference (%)')
avg = mean(y_axis);
txt = ['Mean relative difference: ' num2str(avg,2) '%'];
text(55,max(y_axis),txt,'VerticalAlignment','top')

output = [];