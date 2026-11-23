function weights_layer = computeLayerWeights_HCSMM(k3_layer, U_layer, D_layer, d_layer, fq, cp, method_LinearSolver)
%Compute layer weights of partial waves using HCSMM method.
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

num_layer = length(k3_layer);
xi = 2*pi*fq/cp;

h_layer = cell(1,num_layer);
H_stack = cell(1,num_layer);
for i = 1:num_layer
    k3 = k3_layer{i}; U = U_layer{i}; D = D_layer{i}; d = d_layer(i);
    
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

Uval_Lnbot = solveLinearHomoEqSet(H_stack{end}(4:6,4:6), method_LinearSolver);
Uval_L1top = H_stack{end}(1:3,4:6)*Uval_Lnbot; Sval_L1top = [0; 0; 0];  % sigma33(x3=0)=0, sigma23(x3=0)=0, sigma13(x3=0)=0
USval_L1top = [Uval_L1top; Sval_L1top];

USval_top_L = cell(1,num_layer);
for i = 1:num_layer
    if i==1
        USval_top_L{1} = USval_L1top;
    else
        H1_TT = H_stack{i-1}(1:3,1:3); H1_TB = H_stack{i-1}(1:3,4:6);
        H1_BT = H_stack{i-1}(4:6,1:3); H1_BB = H_stack{i-1}(4:6,4:6);
        
        Uval_bot_prelayer = H1_TB^(-1)*Uval_L1top;    % Uval_bot_L{i-1}
        Sval_bot_prelayer = H1_BB*Uval_bot_prelayer;  % Sval_bot_L{i-1}
        
        USval_top_L{i} = [Uval_bot_prelayer; Sval_bot_prelayer];  % continuity condition: USval_top_L{i} = USval_bot_L{i-1}
    end
end

weights_layer = cell(1,num_layer);
for i = 1:num_layer
    U = U_layer{i};
    D = D_layer{i};
    weights_layer{i} = [U; D]^(-1)*USval_top_L{i};
end

end