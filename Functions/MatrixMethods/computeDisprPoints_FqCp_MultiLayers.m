function [fq_cdd, cp_cdd, fq2cp, cp2fq] = computeDisprPoints_FqCp_MultiLayers(...
    LayerPro, prop_angle, fq_min, fq_max, num_fq, cp_min, cp_max, num_cp, options)
%Compute the dispersion points of multi-layered waveguide using matrix-based methods.
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

% set default parameters for options---------------------------------------
options_default.MatrixModel = 'SMM';  % 4 options: 'TMM', 'SMM', 'HCSMM', 'DVP'
options_default.EQSolver = 'MRCM1D';  % 2 options: 'MRCM1D' and 'PHASE'
options_default.step_cp = (cp_max-cp_min)/num_cp/2;  % unit: m/s
options_default.step_fq = (fq_max-fq_min)/num_fq/2;  % unit: Hz
if nargin<9
    options = options_default;
else
    options = setOptions(options, options_default);
end
MatrixModel = options.MatrixModel;
EQSolver = options.EQSolver;
step_cp = options.step_cp;
step_fq = options.step_fq;
% end to set default parameters for options--------------------------------

num_layer = length(LayerPro);
if num_layer==1
    error('The current structure is a single layer plate, while the multi-layer model is called.');
end

for i = 1:num_layer
    LayerPro(i).prop_angle = prop_angle;  % propagation angle, unit: degree
    LayerPro(i).rot_angle = LayerPro(i).prop_angle - LayerPro(i).orient_angle;  % rotation angle, unit: degree
    LayerPro(i).rotC = StiffnessMatrixRotation(LayerPro(i).C, deg2rad(LayerPro(i).rot_angle));
    LayerPro(i).class_rotC = distinguishStiffnessMatrixClass(LayerPro(i).rotC);
end

FQ_range = linspace(fq_min, fq_max, num_fq);
CP_range = linspace(cp_min, cp_max, num_cp);

disp(['Computing the dispersion points in (fq,cp) domain, using ',MatrixModel,' method...']);

%% loop over fq range, i.e., from fq to cp
disp('Loop over fq...');

parfor_progress(num_fq);
parfor i = 1:num_fq
    fq = FQ_range(i);

    switch MatrixModel
        case 'TMM'
            Fun_cp = @(cp) DispersionEQ_TMM(  LayerPro, fq, cp);
        case 'SMM'
            Fun_cp = @(cp) DispersionEQ_SMM(  LayerPro, fq, cp);
        case 'HCSMM'
            Fun_cp = @(cp) DispersionEQ_HCSMM(LayerPro, fq, cp);
        case 'DVP'
            Fun_cp = @(cp) DispersionEQ_DVP(  LayerPro, fq, cp);
        otherwise
            error('There is something wrong for MatrixModel.');
    end

    switch EQSolver
        case 'MRCM1D'
            cp_roots = Solver_MRCM1D(Fun_cp, cp_min, cp_max, step_cp);
        case 'PHASE'
            cp_roots = Solver_BISECTION_PhaseChange(Fun_cp, cp_min, cp_max, step_cp);
        otherwise
            error('There is something wrong for EQSolver.');
    end

    cp_roots = reshape(cp_roots,1,[]);

    Fq_fq2cp{i} = fq*ones(size(cp_roots));  % cell class
    Cp_fq2cp{i} = cp_roots;                 % cell class

    try parfor_progress; catch, continue; end
end
parfor_progress(0); pause(0.1);

% delete the empty elements in Cp_fq2cp
idx = cellfun('isempty',Cp_fq2cp); Cp_fq2cp(idx) = []; Fq_fq2cp(idx) = [];
% transform cell class into row vector
fq2cp.fq = cell2mat(Fq_fq2cp);
fq2cp.cp = cell2mat(Cp_fq2cp);

%% loop over cp range, i.e., from cp to fq
disp('Loop over cp...');

parfor_progress(num_cp);
parfor i = 1:num_cp
    cp = CP_range(i);

    switch MatrixModel
        case 'TMM'
            Fun_fq = @(fq) DispersionEQ_TMM(LayerPro, fq, cp);
        case 'SMM'
            Fun_fq = @(fq) DispersionEQ_SMM(LayerPro, fq, cp);
        case 'HCSMM'
            Fun_fq = @(fq) DispersionEQ_HCSMM(LayerPro, fq, cp);
        case 'DVP'
            Fun_fq = @(fq) DispersionEQ_DVP(LayerPro, fq, cp);
        otherwise
            error('There is something wrong for MatrixModel.');
    end

    switch EQSolver
        case 'MRCM1D'
            fq_roots = Solver_MRCM1D(Fun_fq, fq_min, fq_max, step_fq);
        case 'PHASE'
            fq_roots = Solver_BISECTION_PhaseChange(Fun_fq, fq_min, fq_max, step_fq);
        otherwise
            error('There is something wrong for EQSolver.');
    end

    fq_roots = reshape(fq_roots,1,[]);

    Cp_cp2fq{i} = cp*ones(size(fq_roots));  % cell class
    Fq_cp2fq{i} = fq_roots;                 % cell class

    try parfor_progress; catch, continue; end
end
parfor_progress(0); pause(0.1);

% delete the empty elements in Fq_cp2fq
idx = cellfun('isempty',Fq_cp2fq); Cp_cp2fq(idx) = []; Fq_cp2fq(idx) = [];
% transform cell class into row vector
cp2fq.fq = cell2mat(Fq_cp2fq);
cp2fq.cp = cell2mat(Cp_cp2fq);

%% concatenate both parts of solutions as the candidate points in the dispersion curves
fq_cdd = [fq2cp.fq, cp2fq.fq];
cp_cdd = [fq2cp.cp, cp2fq.cp];

end