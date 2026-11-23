% Figure 4(b), plot the dispersion curves of OpenGuidedWaves, using DVP method, prop_angle = 0 degree
clc; clear; close all;
addpath(genpath('Functions'));

%% define material properties, data from the caption of Figure 4
C11 = 130; C12 = 6.1; C13 = 6.1; C23 = 5.2; C22 = 11.2; C33 = 11.2; C44 = 3.0; C55 = 4.2; C66 = 4.2;  % unit: Gpa
C = [C11, C12, C13, 0, 0, 0;
    C12, C22, C23, 0, 0, 0;
    C13, C23, C33, 0, 0, 0;
    0, 0, 0, C44, 0, 0;
    0, 0, 0, 0, C55, 0;
    0, 0, 0, 0, 0, C66]*1e9;  % stiffness matrix of orthotropic material, unit: Pa
rho = 1571;   % mass density, unit: kg/m^3
d = 2e-3/16;  % layer thickness, unit: m

%% define layer properties, 16 layers with layup [45/0/-45/90/-45/0/45/90]_s
% layer 1
LayerPro(1).C = C;            % unit: Pa
LayerPro(1).rho = rho;        % unit: kg/m^3
LayerPro(1).d = d;            % unit: m
LayerPro(1).orient_angle =45; % unit: degree

% layer 2
LayerPro(2).C = C;            % unit: Pa
LayerPro(2).rho = rho;        % unit: kg/m^3
LayerPro(2).d = d;            % unit: m
LayerPro(2).orient_angle = 0; % unit: degree

% layer 3
LayerPro(3).C = C;            % unit: Pa
LayerPro(3).rho = rho;        % unit: kg/m^3
LayerPro(3).d = d;            % unit: m
LayerPro(3).orient_angle =-45;% unit: degree

% layer 4
LayerPro(4).C = C;            % unit: Pa
LayerPro(4).rho = rho;        % unit: kg/m^3
LayerPro(4).d = d;            % unit: m
LayerPro(4).orient_angle =90; % unit: degree

% layer 5
LayerPro(5).C = C;            % unit: Pa
LayerPro(5).rho = rho;        % unit: kg/m^3
LayerPro(5).d = d;            % unit: m
LayerPro(5).orient_angle =-45;% unit: degree

% layer 6
LayerPro(6).C = C;            % unit: Pa
LayerPro(6).rho = rho;        % unit: kg/m^3
LayerPro(6).d = d;            % unit: m
LayerPro(6).orient_angle = 0; % unit: degree

% layer 7
LayerPro(7).C = C;            % unit: Pa
LayerPro(7).rho = rho;        % unit: kg/m^3
LayerPro(7).d = d;            % unit: m
LayerPro(7).orient_angle =45; % unit: degree

% layer 8
LayerPro(8).C = C;            % unit: Pa
LayerPro(8).rho = rho;        % unit: kg/m^3
LayerPro(8).d = d;            % unit: m
LayerPro(8).orient_angle =90; % unit: degree

% layer 9
LayerPro(9).C = C;            % unit: Pa
LayerPro(9).rho = rho;        % unit: kg/m^3
LayerPro(9).d = d;            % unit: m
LayerPro(9).orient_angle =90; % unit: degree

% layer 10
LayerPro(10).C = C;            % unit: Pa
LayerPro(10).rho = rho;        % unit: kg/m^3
LayerPro(10).d = d;            % unit: m
LayerPro(10).orient_angle =45; % unit: degree

% layer 11
LayerPro(11).C = C;            % unit: Pa
LayerPro(11).rho = rho;        % unit: kg/m^3
LayerPro(11).d = d;            % unit: m
LayerPro(11).orient_angle = 0; % unit: degree

% layer 12
LayerPro(12).C = C;            % unit: Pa
LayerPro(12).rho = rho;        % unit: kg/m^3
LayerPro(12).d = d;            % unit: m
LayerPro(12).orient_angle =-45;% unit: degree

% layer 13
LayerPro(13).C = C;            % unit: Pa
LayerPro(13).rho = rho;        % unit: kg/m^3
LayerPro(13).d = d;            % unit: m
LayerPro(13).orient_angle =90; % unit: degree

% layer 14
LayerPro(14).C = C;            % unit: Pa
LayerPro(14).rho = rho;        % unit: kg/m^3
LayerPro(14).d = d;            % unit: m
LayerPro(14).orient_angle =-45;% unit: degree

% layer 15
LayerPro(15).C = C;            % unit: Pa
LayerPro(15).rho = rho;        % unit: kg/m^3
LayerPro(15).d = d;            % unit: m
LayerPro(15).orient_angle = 0; % unit: degree

% layer 16
LayerPro(16).C = C;            % unit: Pa
LayerPro(16).rho = rho;        % unit: kg/m^3
LayerPro(16).d = d;            % unit: m
LayerPro(16).orient_angle =45; % unit: degree

num_layer = length(LayerPro);

TH = 0;
for i = 1:num_layer
    TH = TH + LayerPro(i).d;
end

%% define solving parameters
prop_angle = 0;  % propagation angle, unit: degree

fq_min = 1e3; fq_max = 2e6; num_fq = 300;  % unit: Hz
step_fq = (fq_max-fq_min)/num_fq;

cp_min = 10; cp_max = 10e3; num_cp = 300;  % unit: m/s
step_cp = (cp_max-cp_min)/num_cp;

%% solve
tic
[fq_cdd, cp_cdd] = computeDisprPoints_FqCp_MultiLayers(...
    LayerPro, prop_angle, fq_min, fq_max, num_fq, cp_min, cp_max, num_cp,...
    struct('MatrixModel','DVP'));
toc

figure; set(gcf, 'windowstyle','docked');
plot(fq_cdd/1e3, cp_cdd/1e3, 'k.'); box on; grid on;
xlabel('Frequency (kHz)'); ylabel('Phase velocity (km/s)');
text(80, 9.2, '(b) DVP', 'fontsize',12);

%% dense the veering regions
% region_dense(1,:) = [fq_min, fq_max, num_fq, step_fq, cp_min, cp_max, num_cp, step_cp];
region_dense = [
    0.9e3, 8e3, 10, 1e3, 140, 450, 10, 10;         % region 1
    2.2e5, 2.4e5, 10, 1e3, 1280, 1310, 10, 1;      % region 2
    3.2e5, 3.4e5, 10, 1e3, 1325, 1340, 10, 1;      % region 3
    4.0e5, 4.2e5, 10, 1e3, 1345, 1355, 10, 1;      % region 4
    4.4e5, 4.6e5, 10, 1e3, 1355, 1362, 10, 0.5;    % region 5
    5.85e5, 6.05e5, 10, 1e3, 1372, 1378, 10, 0.1;  % region 6
    6.40e5, 6.60e5, 10, 1e3, 1375, 1385, 10, 0.1;  % region 7
    7.60e5, 7.80e5, 10, 1e3, 1385, 1395, 10, 0.1;  % region 8
    8.05e5, 8.25e5, 10, 1e3, 1392, 1394, 10, 0.1;  % region 9
    8.35e5, 8.65e5, 10, 1e3, 1394, 1397, 10, 0.1;  % region 10
    9.35e5, 9.55e5, 10, 1e3, 1400, 1403, 10, 0.1;  % region 11
    1.00e6, 1.02e6, 10, 1e3, 1404, 1406, 10, 0.1;  % region 12
    1.06e6, 1.08e6, 10, 1e3, 1407, 1410, 10, 0.1;  % region 13
    1.15e6, 1.19e6, 20, 1e3, 1412, 1416, 10, 0.1;  % region 14
    1.25e6, 1.27e6, 10, 1e3, 1416, 1419, 10, 0.1;  % region 15
    1.30e6, 1.39e6, 20, 1e3, 1419, 1425, 10, 0.05; % region 16
    1.52e6, 1.56e6, 20, 1e3, 1428, 1431, 10, 0.1;  % region 17
    1.615e6, 1.635e6, 10, 1e3, 1432, 1433, 10, 0.05;     % region 18
    1.665e6, 1.680e6, 10, 1e3, 1434.0, 1434.5, 10, 0.01; % region 19
    1.74e6, 1.76e6, 10, 1e3, 1436, 1438, 10, 0.05;       % region 20
    1.86e6, 1.88e6, 10, 1e3, 1440, 1442, 10, 0.05;       % region 21
    1.965e6, 1.999e6, 20, 1e3, 1443, 1445, 10, 0.05;     % region 22
    ...
    9.20e5, 9.65e5, 20, 1e3, 1580, 1625, 10, 0.1;    % region 23
    1.00e6, 1.13e6, 20, 1e3, 1514, 1569, 10, 0.1;    % region 24
    1.16e6, 1.27e6, 20, 1e3, 1488, 1511, 10, 0.1;    % region 25
    1.31e6, 1.47e6, 10, 1e3, 1469, 1484, 10, 0.1;    % region 26
    1.500e6, 1.535e6, 20, 1e3, 1466, 1469, 10, 0.05; % region 27
    1.59e6, 1.69e6, 20, 1e3, 1462, 1465, 10, 0.02;   % region 28
    1.75e6, 1.77e6, 10, 1e3, 1461, 1462, 10, 0.01;   % region 29
    1.935e6, 1.955e6, 10, 1e3, 1460, 1461, 10, 0.01; % region 30
    1.97e6, 1.999e6, 40, 1e3, 1619, 1628, 10, 0.05;  % region 31
    ...
    6.20e5, 6.60e5, 30, 1e3, 3000, 3300, 30, 1;      % region 32
    ];
figure; hold on;grid on;box on; set(gcf, 'windowstyle','docked');
plot(fq_cdd, cp_cdd, 'k.');
for i = 1:size(region_dense,1)
    fq_min_dense = region_dense(i,1); fq_max_dense = region_dense(i,2); cp_min_dense = region_dense(i,5); cp_max_dense = region_dense(i,6);
    plot([fq_min_dense,fq_max_dense,fq_max_dense,fq_min_dense,fq_min_dense], [cp_min_dense,cp_min_dense,cp_max_dense,cp_max_dense,cp_min_dense], 'r-');
    text(fq_min_dense, cp_min_dense, num2str(i), 'color','r');
end
xlabel('Frequency (Hz)'); ylabel('Phase velocity (m/s)');
xlim([0,fq_max]);ylim([0,cp_max]);


clear options
tic
cp_dense = []; fq_dense = [];
for i = 1:size(region_dense,1)
    disp(['Computing region ',num2str(i),'/',num2str(size(region_dense,1)),'...']);
    fq_min_dense = region_dense(i,1); fq_max_dense = region_dense(i,2); num_fq_dense = region_dense(i,3); step_fq_dense = region_dense(i,4);
    cp_min_dense = region_dense(i,5); cp_max_dense = region_dense(i,6); num_cp_dense = region_dense(i,7); step_cp_dense = region_dense(i,8);
    options.step_cp = step_cp_dense; options.step_fq = step_fq_dense; options.MatrixModel = 'DVP';

    [fq_dense_temp, cp_dense_temp] = computeDisprPoints_FqCp_MultiLayers(LayerPro, prop_angle,...
        fq_min_dense, fq_max_dense, num_fq_dense, cp_min_dense, cp_max_dense, num_cp_dense, options);

    cp_dense = [cp_dense, cp_dense_temp]; fq_dense = [fq_dense, fq_dense_temp];
end
toc


figure; hold on;grid on;box on; set(gcf, 'windowstyle','docked');
plot(fq_cdd, cp_cdd, 'k.');
for i = 1:size(region_dense,1)
    fq_min_dense = region_dense(i,1); fq_max_dense = region_dense(i,2); cp_min_dense = region_dense(i,5); cp_max_dense = region_dense(i,6);
    plot([fq_min_dense,fq_max_dense,fq_max_dense,fq_min_dense,fq_min_dense], [cp_min_dense,cp_min_dense,cp_max_dense,cp_max_dense,cp_min_dense], 'r-');
    text(fq_min_dense, cp_min_dense, num2str(i), 'color','r');
end
plot(fq_dense, cp_dense, 'r.');
xlabel('Frequency (Hz)'); ylabel('Phase velocity (m/s)');
xlim([0,fq_max]);ylim([0,cp_max]);

% concatenate both parts of solutions
fq_cdd = [fq_cdd, fq_dense]; cp_cdd = [cp_cdd, cp_dense];

%% plot Figure 4(b) DVP
fig4b = figure;
plot(fq_cdd/1e3, cp_cdd/1e3, 'k.'); box on;
xlabel('Frequency (kHz)'); ylabel('Phase velocity (km/s)');
set(gca, 'xtick',[0,250,500,750,1000,1250,1500,1750,2000], 'ytick',[0,2,4,6,8,10], 'fontsize',12);
set(gca, 'position',[0.085, 0.12, 0.88, 0.86]);
text(80, 9.2, '(b) DVP', 'fontsize',12);

% save this figure
% exportgraphics(fig4b, 'Data_Images/Figure_4b.pdf', 'ContentType','vector');

%% save data for direct later usage
% save('Data_Images/data_Fig4b_OpenGuidedWaves_FqCp_0deg_DVP.mat','fq_cdd','cp_cdd');
% load('Data_Images/data_Fig4b_OpenGuidedWaves_FqCp_0deg_DVP.mat');
