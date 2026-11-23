function [y, h_layer, H_stack] = DispersionEQ_HCSMM(LayerPro, fq, cp)
%Construct the dispersion equation using hybrid compliance-stiffness matrix method (HCSMM).
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

num_layer = length(LayerPro);
% k1 = 2*pi*fq/cp;  % wavenumber in x1 direction

h_layer = cell(1,num_layer);
H_stack = cell(1,num_layer);
for i = 1:num_layer
    d = LayerPro(i).d;
    [k3, U, D] = solveStrohMatrix(LayerPro(i).rotC, LayerPro(i).rho, fq, cp);
    
    k3_T = k3(1:3);                   k3_B = k3(4:6);
    Lambda_T = diag(exp(+1i*k3_T*d)); Lambda_B = diag(exp(-1i*k3_B*d));
    U_T = U(:,1:3);                   U_B = U(:,4:6);
    D_T = D(:,1:3);                   D_B = D(:,4:6);
    
    h_temp1 = [U_T, U_B*Lambda_B; D_T*Lambda_T, D_B];
    h_temp2 = [D_T, D_B*Lambda_B; U_T*Lambda_T, U_B];
    h_layer{i} = h_temp1*h_temp2^(-1);
    
    % the recursive process
    if i==1
        H_stack{1} = h_layer{1};
    else
        h_TT = h_layer{i}(1:3,1:3); h_TB = h_layer{i}(1:3,4:6);
        h_BT = h_layer{i}(4:6,1:3); h_BB = h_layer{i}(4:6,4:6);
        
        H1_TT = H_stack{i-1}(1:3,1:3); H1_TB = H_stack{i-1}(1:3,4:6);
        H1_BT = H_stack{i-1}(4:6,1:3); H1_BB = H_stack{i-1}(4:6,4:6);
        
        delta_inv = (eye(3)-h_TT*H1_BB)^(-1);
        
        H_stack{i} = [H1_TT+H1_TB*delta_inv*h_TT*H1_BT, H1_TB*delta_inv*h_TB;
            h_BT*(eye(3)+H1_BB*delta_inv*h_TT)*H1_BT, h_BT*H1_BB*delta_inv*h_TB+h_BB];
    end
end

y = det( H_stack{end}(4:6,4:6) );

end