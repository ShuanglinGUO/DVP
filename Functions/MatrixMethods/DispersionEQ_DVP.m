function [y, s_layer, S_stack] = DispersionEQ_DVP(LayerPro, fq, cp)
%Construct the dispersion equation using dual-variable and position method (DVP).
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

num_layer = length(LayerPro);
% k1 = 2*pi*fq/cp;  % wavenumber in x1 direction

s_layer = cell(1,num_layer);
S_stack = cell(1,num_layer);
for i = 1:num_layer
    d = LayerPro(i).d;
    [k3, U, D] = solveStrohMatrix(LayerPro(i).rotC, LayerPro(i).rho, fq, cp);
    
    k3_T = k3(1:3);                   k3_B = k3(4:6);
    Lambda_T = diag(exp(+1i*k3_T*d)); Lambda_B = diag(exp(-1i*k3_B*d));
    U_T = U(:,1:3);                   U_B = U(:,4:6);
    D_T = D(:,1:3);                   D_B = D(:,4:6);
    
    s_temp1 = [U_T, U_B*Lambda_B; D_T*Lambda_T, D_B];
    s_temp2 = [U_T*Lambda_T, U_B; D_T, D_B*Lambda_B];
    s_layer{i} = s_temp1*s_temp2^(-1);
    
    % the recursive process
    if i==1
        S_stack{1} = s_layer{1};
    else
        s_TB = s_layer{i}(1:3,1:3); s_TT = s_layer{i}(1:3,4:6);
        s_BB = s_layer{i}(4:6,1:3); s_BT = s_layer{i}(4:6,4:6);
        
        S1_TB = S_stack{i-1}(1:3,1:3); S1_TT = S_stack{i-1}(1:3,4:6);
        S1_BB = S_stack{i-1}(4:6,1:3); S1_BT = S_stack{i-1}(4:6,4:6);
        
        delta_inv = (eye(3)-s_TT*S1_BB)^(-1);
        
        S_stack{i} = [S1_TB*delta_inv*s_TB, S1_TT+S1_TB*delta_inv*s_TT*S1_BT;
            s_BB+s_BT*S1_BB*delta_inv*s_TB, s_BT*(eye(3)+S1_BB*delta_inv*s_TT)*S1_BT];
    end
end

y = det( S_stack{end}(4:6,1:3) );

end