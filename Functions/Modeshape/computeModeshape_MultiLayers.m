function [x3_vec, modeshape_disp_cdd, modeshape_stress_cdd] = computeModeshape_MultiLayers(LayerPro, prop_angle, fq_cdd, cp_cdd, options)
%Compute (displacement and stress) modeshapes for all (fq_cdd,cp_cdd) points.
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

% set default parameters for options---------------------------------------
options_default.num_x3 = 100;
options_default.Complex2RealMethod = 'Rotation';  % 3 options: 'None', 'Rotation' and 'RealImag'
options_default.MatrixModel = 'SMM';              % 4 options: 'TMM', 'SMM', 'HCSMM', 'DVP'
options_default.SwitchNorm  = 'ON';               % 2 options: 'ON' and 'OFF'
options_default.method_LinearSolver = 'eig';      % 4 options: 'indeterminacy', 'normalization', 'eig', 'svd' ('normalization' and 'eig' are two preferred methods)
if nargin<5
    options = options_default;
else
    options = setOptions(options, options_default);
end
num_x3 = options.num_x3;
Complex2RealMethod = options.Complex2RealMethod;
MatrixModel = options.MatrixModel;
SwitchNorm  = options.SwitchNorm;
method_LinearSolver = options.method_LinearSolver;
% end to set default parameters for options--------------------------------


num_layer  = length(LayerPro);
num_points = length(fq_cdd);

for i = 1:num_layer
    LayerPro(i).prop_angle = prop_angle;  % propagation angle, unit: degree
    LayerPro(i).rot_angle = LayerPro(i).prop_angle - LayerPro(i).orient_angle;  % rotation angle, unit: degree
    LayerPro(i).rotC = StiffnessMatrixRotation(LayerPro(i).C, deg2rad(LayerPro(i).rot_angle));
    LayerPro(i).class_rotC = distinguishStiffnessMatrixClass(LayerPro(i).rotC);
end

modeshape_disp_cdd   = zeros(3, num_layer*num_x3, num_points);
modeshape_stress_cdd = zeros(3, num_layer*num_x3, num_points);
for n = 1:num_points    % this for-loop is enable for computing all (fq_cdd, cp_cdd) points
    fq = fq_cdd(n); cp = cp_cdd(n);
    
    k3_layer = cell(1,num_layer);
    U_layer = cell(1,num_layer);
    D_layer = cell(1,num_layer);
    d_layer = zeros(1,num_layer);
    for i = 1:num_layer
        [k3_layer{i}, U_layer{i}, D_layer{i}] = solveStrohMatrix(LayerPro(i).rotC, LayerPro(i).rho, fq, cp);
        d_layer(i) = LayerPro(i).d;
    end
    
    switch MatrixModel
        case 'TMM'
            weights_layer = computeLayerWeights_TMM(  k3_layer, U_layer, D_layer, d_layer, fq, cp, method_LinearSolver);
        case 'SMM'
            weights_layer = computeLayerWeights_SMM(  k3_layer, U_layer, D_layer, d_layer, fq, cp, method_LinearSolver);
        case 'HCSMM'
            weights_layer = computeLayerWeights_HCSMM(k3_layer, U_layer, D_layer, d_layer, fq, cp, method_LinearSolver);
        case 'DVP'
            weights_layer = computeLayerWeights_DVP(  k3_layer, U_layer, D_layer, d_layer, fq, cp, method_LinearSolver);
        otherwise
            error('Invalid value of options.MatrixModel.');
    end
    
    d_cumsum = 0;
    x3_layer = cell(1,num_layer);
    modeshape_disp_layer   = cell(1,num_layer);
    modeshape_stress_layer = cell(1,num_layer);
    for i = 1:num_layer
        d = LayerPro(i).d;
        k3 = k3_layer{i}; U = U_layer{i}; D = D_layer{i}; weights = weights_layer{i};
        
        x3 = linspace(0, d, num_x3);  % local coordinate of thickness, unit: m
        for j = 1:num_x3
            diagE = diag(exp( 1i*k3*x3(j) ));
            modeshape_disp_layer{i}(:,j)   = U*diagE*weights;
            modeshape_stress_layer{i}(:,j) = D*diagE*weights;
        end
        x3_layer{i} = x3 + d_cumsum; d_cumsum = d_cumsum + d;  % cumsum coordinate of thickness, unit: m
    end
    
    modeshape_disp   = cell2mat(modeshape_disp_layer);
    modeshape_stress = cell2mat(modeshape_stress_layer);
    x3_vec = cell2mat(x3_layer);
    
    switch Complex2RealMethod
        case 'None'
            %modeshape_disp   = modeshape_disp;
            %modeshape_stress = modeshape_stress;
        case 'Rotation'
            modeshape_disp   = modeshape_complex2real_RotationMethod(modeshape_disp);
            modeshape_stress = modeshape_complex2real_RotationMethod(modeshape_stress);
        case 'RealImag'
            modeshape_disp   = real(modeshape_disp)   + imag(modeshape_disp);
            modeshape_stress = real(modeshape_stress) + imag(modeshape_stress);
        otherwise
            error('Complex2RealMethod has only three cases: ''None'', ''Rotation'', ''RealImag''.');
    end
    
    if strcmp(SwitchNorm, 'ON')
        modeshape_disp   = modeshape_normalization(modeshape_disp);
        modeshape_stress = modeshape_normalization(modeshape_stress);
    end
    
    modeshape_disp_cdd(:,:,n)   = modeshape_disp;
    modeshape_stress_cdd(:,:,n) = modeshape_stress;
end
end

function modeshape_real = modeshape_complex2real_RotationMethod(modeshape_complex)
[num_disp, num_x3] = size(modeshape_complex);
modeshape_complex_concatenated = reshape(modeshape_complex.', 1, []);

% rotate modeshape to make the second complex number to be real, taking the second point may lead to stable modeshape
modeshape_rotated_temp = modeshape_complex_concatenated.*exp( -1j*angle(modeshape_complex_concatenated(2)) );

modeshape_rotated = zeros(num_disp, num_x3);
for i = 1:num_disp
    modeshape_temp = modeshape_rotated_temp( (i-1)*num_x3+1:i*num_x3 );
    modeshape_rotated(i,:) = modeshape_dominant( modeshape_temp );
end
modeshape_real = modeshape_rotated;
end

function modeshape = modeshape_dominant(modeshape_complex)
% judge which part (real or imag) is dominant in a complex-valued modeshape
if max(abs(real(modeshape_complex))) >= max(abs(imag(modeshape_complex)))
    modeshape = real(modeshape_complex);
elseif max(abs(real(modeshape_complex))) < max(abs(imag(modeshape_complex)))
    modeshape = imag(modeshape_complex);
else
    error('There is something wrong!');
end
end

function modeshape_norm = modeshape_normalization(modeshape)
[max_val, idx_max] = max(abs(modeshape(:)));
sign_max = sign( modeshape(idx_max) );
modeshape_norm = modeshape/max_val*sign_max;  % normalized modeshape
end