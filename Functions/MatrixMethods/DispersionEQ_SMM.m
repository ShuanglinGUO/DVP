function [y, k_layer, K_stack] = DispersionEQ_SMM(LayerPro, fq, cp)
%Construct the dispersion equation using stiffness matrix method (SMM).
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

num_layer = length(LayerPro);
% k1 = 2*pi*fq/cp;  % wavenumber in x1 direction

k_layer = cell(1,num_layer);
K_stack = cell(1,num_layer);
for i = 1:num_layer
    d = LayerPro(i).d;
    [k3, U, D] = solveStrohMatrix(LayerPro(i).rotC, LayerPro(i).rho, fq, cp);
    
    k3_T = k3(1:3);                   k3_B = k3(4:6);
    Lambda_T = diag(exp(+1i*k3_T*d)); Lambda_B = diag(exp(-1i*k3_B*d));
    U_T = U(:,1:3);                   U_B = U(:,4:6);
    D_T = D(:,1:3);                   D_B = D(:,4:6);
    
    kappa_u = [U_T, U_B*Lambda_B; U_T*Lambda_T, U_B];
    kappa_sigma = [D_T, D_B*Lambda_B; D_T*Lambda_T, D_B];
    k_layer{i} = kappa_sigma*kappa_u^(-1);
    
    % the recursive process
    if i==1
        K_stack{1} = k_layer{1};
    else
        k_TT = k_layer{i}(1:3,1:3); k_TB = k_layer{i}(1:3,4:6);
        k_BT = k_layer{i}(4:6,1:3); k_BB = k_layer{i}(4:6,4:6);
        
        K1_TT = K_stack{i-1}(1:3,1:3); K1_TB = K_stack{i-1}(1:3,4:6);
        K1_BT = K_stack{i-1}(4:6,1:3); K1_BB = K_stack{i-1}(4:6,4:6);
        
        mat_inv = (k_TT-K1_BB)^(-1);
        
        K_stack{i} = [K1_TT+K1_TB*mat_inv*K1_BT, -K1_TB*mat_inv*k_TB;
            k_BT*mat_inv*K1_BT, k_BB-k_BT*mat_inv*k_TB];
    end
end

y = det( K_stack{end} );

end