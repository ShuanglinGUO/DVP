% Figure 4(a), plot the dispersion curves of OpenGuidedWaves, using TMM method, prop_angle = 0 degree
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

fq_min = 1e3; fq_max = 2e6; num_fq = 200;  % unit: Hz
step_fq = (fq_max-fq_min)/num_fq;

cp_min = 10; cp_max = 10e3; num_cp = 200;  % unit: m/s
step_cp = (cp_max-cp_min)/num_cp;

%% solve
tic
[fq_cdd, cp_cdd] = computeDisprPoints_FqCp_MultiLayers(...
    LayerPro, prop_angle, fq_min, fq_max, num_fq, cp_min, cp_max, num_cp,...
    struct('MatrixModel','TMM'));
toc

%% plot Figure 4(a) TMM
fig4a = figure;
plot(fq_cdd/1e3, cp_cdd/1e3, 'k.'); box on;
xlabel('Frequency (kHz)'); ylabel('Phase velocity (km/s)');
set(gca, 'xtick',[0,250,500,750,1000,1250,1500,1750,2000], 'ytick',[0,2,4,6,8,10], 'fontsize',12);
set(gca, 'position',[0.085, 0.12, 0.88, 0.86]);
text(80, 9.2, '(a) TMM', 'fontsize',12);

% save this figure
% exportgraphics(fig4a, 'Data_Images/Figure_4a.pdf', 'ContentType','vector');

%% save data for direct later usage
% save('Data_Images/data_Fig4a_OpenGuidedWaves_FqCp_0deg_TMM.mat','fq_cdd','cp_cdd');
% load('Data_Images/data_Fig4a_OpenGuidedWaves_FqCp_0deg_TMM.mat');
