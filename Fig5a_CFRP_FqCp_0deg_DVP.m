% Figure 5(a), plot the dispersion curves of CFRP laminate, using DVP method, prop_angle = 0 degree
clc; clear; close all;
addpath(genpath('Functions'));

%% define material properties, data from the caption of Figure 5
E1 = 143e9; E2 = 9.1e9; G12 = 4.8e9; nu12 = 0.3; nu23 = 0.3;
rho = 1650;    % unit: kg/m^3
d = 0.17e-3;   % unit: m

G23 = E2/2/(1+nu23);

% The stiffness matrix in the principal direction
C = [   1/E1, -nu12/E1, -nu12/E1,     0,     0,     0;
    -nu12/E1,     1/E2, -nu23/E2,     0,     0,     0;
    -nu12/E1, -nu23/E2,     1/E2,     0,     0,     0;
           0,        0,        0, 1/G23,     0,     0;
           0,        0,        0,     0, 1/G12,     0;
           0,        0,        0,     0,     0, 1/G12]^(-1);  % unit: Pa

%% define layer properties, 8 layers with layup [0/90/+45/-45]_s
% layer 1
LayerPro(1).C = C;            % unit: Pa
LayerPro(1).rho = rho;        % unit: kg/m^3
LayerPro(1).d = d;            % unit: m
LayerPro(1).orient_angle = 0; % unit: degree

% layer 2
LayerPro(2).C = C;            % unit: Pa
LayerPro(2).rho = rho;        % unit: kg/m^3
LayerPro(2).d = d;            % unit: m
LayerPro(2).orient_angle = 90;% unit: degree

% layer 3
LayerPro(3).C = C;            % unit: Pa
LayerPro(3).rho = rho;        % unit: kg/m^3
LayerPro(3).d = d;            % unit: m
LayerPro(3).orient_angle = 45;% unit: degree

% layer 4
LayerPro(4).C = C;            % unit: Pa
LayerPro(4).rho = rho;        % unit: kg/m^3
LayerPro(4).d = d;            % unit: m
LayerPro(4).orient_angle =-45;% unit: degree

% layer 5
LayerPro(5).C = C;            % unit: Pa
LayerPro(5).rho = rho;        % unit: kg/m^3
LayerPro(5).d = d;            % unit: m
LayerPro(5).orient_angle =-45;% unit: degree

% layer 6
LayerPro(6).C = C;            % unit: Pa
LayerPro(6).rho = rho;        % unit: kg/m^3
LayerPro(6).d = d;            % unit: m
LayerPro(6).orient_angle = 45;% unit: degree

% layer 7
LayerPro(7).C = C;            % unit: Pa
LayerPro(7).rho = rho;        % unit: kg/m^3
LayerPro(7).d = d;            % unit: m
LayerPro(7).orient_angle = 90;% unit: degree

% layer 8
LayerPro(8).C = C;            % unit: Pa
LayerPro(8).rho = rho;        % unit: kg/m^3
LayerPro(8).d = d;            % unit: m
LayerPro(8).orient_angle = 0; % unit: degree

num_layer = length(LayerPro);

TH = 0;
for i = 1:num_layer
    TH = TH + LayerPro(i).d;
end

%% define solving parameters
prop_angle = 0;  % propagation angle, unit: degree

fd_min = 1; fd_max = 1000;  % unit: Hz.m
fq_min = fd_min/TH; fq_max = fd_max/TH; num_fq = 100;  % unit: Hz
step_fq = (fq_max-fq_min)/num_fq;

cp_min = 10; cp_max = 7000; num_cp = 100;  % unit: m/s
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

% save data for direct later usage
% save('Data_Images/data_Fig5a_CFRP_FqCp_0deg_DVP.mat','fq_cdd','cp_cdd');
% load('Data_Images/data_Fig5a_CFRP_FqCp_0deg_DVP.mat');

%% compute modeshape and try to classify symmetry for candidate points (fq_cdd, cp_cdd)
[x3_vec, modeshape_disp_cdd, modeshape_stress_cdd] = computeModeshape_MultiLayers(LayerPro, prop_angle, fq_cdd, cp_cdd,...
    struct('MatrixModel','DVP', 'num_x3',50));

[fq_sym, cp_sym, modeshape_sym,...
    fq_anti, cp_anti, modeshape_anti,...
    fq_none, cp_none, modeshape_none] = classifySymmetry2D(...
    fq_cdd, cp_cdd, modeshape_disp_cdd, struct('tol',1e-1));

figure; set(gcf, 'windowstyle','docked'); hold on;grid on;box on;
plot(fq_none/1e3, cp_none/1e3, 'k.');
plot(fq_sym/1e3,  cp_sym/1e3,  'b.');
plot(fq_anti/1e3, cp_anti/1e3, 'r.');
xlabel('Frequency (kHz)'); ylabel('Phase velocity (km/s)');
xlim([0,fq_max/1e3]); ylim([0,cp_max/1e3]);
title('Sym modes (blue), Anti modes (red), None symmetry modes (black)');

%% extract continuous curves of each mode, and compute group velocity

% extract continuous dispersion curves of S0 mode
idx_S0 = cp_sym>=5000 & cp_sym<=7000;
fq_mode{1} = fq_sym(idx_S0); cp_mode{1} = cp_sym(idx_S0);

% extract continuous dispersion curves of SH0 mode
idx_SH0 = cp_sym>=2500 & cp_sym<=4000;
fq_mode{2} = fq_sym(idx_SH0); cp_mode{2} = cp_sym(idx_SH0);

% extract continuous dispersion curves of A0 mode
idx_A0 = cp_anti<=2000;
fq_mode{3} = fq_anti(idx_A0); cp_mode{3} = cp_anti(idx_A0);

% extract continuous dispersion curves of SH1 mode
idx_SH1 = cp_anti>=3800;
fq_mode{4} = fq_anti(idx_SH1); cp_mode{4} = cp_anti(idx_SH1);

% re-sort each mode based on the increasing frequency
num_mode = length(fq_mode);
for r = 1:num_mode
    [~,idx] = sort(fq_mode{r});
    fq_mode{r} = fq_mode{r}(idx); cp_mode{r} = cp_mode{r}(idx);
end

% plot continuous dispersion curves for each mode
figure; set(gcf, 'windowstyle','docked'); hold on;grid on;box on;
for r = 1:num_mode
    plot( fq_mode{r}/1e3, cp_mode{r}/1e3, 'b-', 'linewidth',2 );
end
xlabel('Frequency (kHz)'); ylabel('Phase velocity (km/s)');
xlim([0,fq_max/1e3]); ylim([0,cp_max/1e3]);

% compute group velocity of each mode, using isotropic model
fq4cg_mode = cell(1,num_mode); cg_mode = cell(1,num_mode);
for r = 1:num_mode
    [fq4cg_mode{r}, cg_mode{r}] = computeGroupVelocity_Isotropic(fq_mode{r}, cp_mode{r});
end

%% plot Figure 5(a), DVP method for CFRP laminate
fig5a = figure; box on;hold on;

hcp = plot(fq_mode{1}*TH/1e3, cp_mode{1}/1e3, 'b-', 'linewidth',2);    % S0 mode, cp
hcp = plot(fq_mode{2}*TH/1e3, cp_mode{2}/1e3, 'b-', 'linewidth',2);    % SH0 mode, cp
hcp = plot(fq_mode{3}*TH/1e3, cp_mode{3}/1e3, 'b-', 'linewidth',2);    % A0 mode, cp
hcp = plot(fq_mode{4}*TH/1e3, cp_mode{4}/1e3, 'b-', 'linewidth',2);    % SH1 mode, cp

hcg = plot(fq4cg_mode{1}*TH/1e3, cg_mode{1}/1e3, 'g:', 'linewidth',2); % S0 mode, cg
hcg = plot(fq4cg_mode{2}*TH/1e3, cg_mode{2}/1e3, 'g:', 'linewidth',2); % SH0 mode, cg
hcg = plot(fq4cg_mode{3}*TH/1e3, cg_mode{3}/1e3, 'g:', 'linewidth',2); % A0 mode, cg
hcg = plot(fq4cg_mode{4}*TH/1e3, cg_mode{4}/1e3, 'g:', 'linewidth',2); % SH1 mode, cg

xlabel('Frequency$\times$Thickness (MHz.mm)', 'interpreter','latex'); ylabel('Wave velocity (km/s)');
xlim([0,0.8]); ylim([0,7]);
set(gca, 'position',[0.069, 0.115, 0.905, 0.865]);
set(gca, 'xtick',[0,0.2,0.4,0.6,0.8,1.0], 'fontsize',12);

% compare with experimental data
load('EXP_data/EXP_Carbon_0degree.mat');
S0 = EXP_Carbon_0degree.S0;
A0 = EXP_Carbon_0degree.A0;
hexp = plot(S0(:,1), S0(:,2), 'ro', 'linewidth',1.5);
hexp = plot(A0(:,1), A0(:,2), 'ro', 'linewidth',1.5);

legend1 = legend([hcp, hcg, hexp], 'Phase velocity', 'Group velocity', 'Exp. data points');
set(legend1, 'Position',[0.6800, 0.1330, 0.2768, 0.1434]);

text(0.05, 6.6, '(a) CFRP Laminate', 'fontsize',12);
text(0.05, 5.8, 'S0 mode', 'fontsize',12);
text(0.05, 3.8, 'SH0 mode', 'fontsize',12);
text(0.05, 1.9, 'A0 mode', 'fontsize',12);
text(0.87, 6.6, 'SH1 mode', 'fontsize',12);

% save this figure
% exportgraphics(fig5a, 'Data_Images/Figure_5a.pdf', 'ContentType','vector');

%% plot modeshape at a picked point
n = 100;  % the picked index
fq_pick = fq_cdd(n);
cp_pick = cp_cdd(n);
xi_pick = 2*pi*fq_pick/cp_pick;

% plot dispersion points and highlight the picked point
hf1 = figure; hold on;grid on;box on; set(hf1,'windowstyle','docked');
plot(fq_cdd/1e3,  cp_cdd/1e3,  'k.');
plot(fq_pick/1e3, cp_pick/1e3, 'bo', 'markerfacecolor','b'); hold off;
xlabel('Frequency (kHz)'); ylabel('Phase velocity (km/s)');

% compute displacement and stress modeshapes at the picked point
[x3_pick, modeshape_disp_pick, modeshape_stress_pick] = computeModeshape_MultiLayers(LayerPro, prop_angle, fq_pick, cp_pick,...
    struct('MatrixModel','DVP', 'num_x3',10));

% plot displacement modeshape
hf2 = figure; hold on;grid on;box on; set(hf2,'windowstyle','docked');
plot( modeshape_disp_pick(1,:), x3_pick, 'k-', 'linewidth',2 );    % u1
plot( modeshape_disp_pick(2,:), x3_pick, 'k:', 'linewidth',2 );    % u2
plot( modeshape_disp_pick(3,:), x3_pick, 'k-.', 'linewidth',2 );   % u3
xlabel('Relative amplitude'); ylabel('Thickness (m)'); legend('u1','u2','u3');
ylim([min(x3_pick),max(x3_pick)]);
title(['Displacement modeshape at point (', num2str(fq_pick/1e3),' kHz, ',num2str(cp_pick/1e3),' km/s)']);

% plot stress modeshape
hf3 = figure; hold on;grid on;box on; set(hf3,'windowstyle','docked');
plot( modeshape_stress_pick(1,:), x3_pick, 'k-', 'linewidth',2 );    % sigma33
plot( modeshape_stress_pick(2,:), x3_pick, 'k:', 'linewidth',2 );    % sigma23
plot( modeshape_stress_pick(3,:), x3_pick, 'k-.', 'linewidth',2 );   % sigma13
xlabel('Relative amplitude'); ylabel('Thickness (m)'); legend('sigma33','sigma23','sigma13');
ylim([min(x3_pick),max(x3_pick)]);
title(['Stress modeshape at point (', num2str(fq_pick/1e3),' kHz, ', num2str(cp_pick/1e3),' km/s)']);

% make modeshape animation
ModeshapeAnimation(x3_pick, modeshape_disp_pick, cp_pick, xi_pick,...
    struct('num_wavelength',2, 'num_period',2, 'scale',1));
