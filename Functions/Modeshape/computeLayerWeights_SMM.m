function weights_layer = computeLayerWeights_SMM(k3_layer, U_layer, D_layer, d_layer, fq, cp, method_LinearSolver)
%Compute layer weights of partial waves using SMM method.
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

num_layer = length(k3_layer);
xi = 2*pi*fq/cp;

k_layer = cell(1,num_layer);
K_stack = cell(1,num_layer);
for i = 1:num_layer
    k3 = k3_layer{i}; U = U_layer{i}; D = D_layer{i}; d = d_layer(i);
    
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

Uval_L1top_Lnbot = solveLinearHomoEqSet(K_stack{end}, method_LinearSolver);
Uval_L1top = Uval_L1top_Lnbot(1:3); Sval_L1top = [0; 0; 0];  % sigma33(x3=0)=0, sigma23(x3=0)=0, sigma13(x3=0)=0
USval_L1top = [Uval_L1top; Sval_L1top];

USval_top_L = cell(1,num_layer);
for i = 1:num_layer
    if i==1
        USval_top_L{1} = USval_L1top;
    else
        K1_TT = K_stack{i-1}(1:3,1:3); K1_TB = K_stack{i-1}(1:3,4:6); 
        K1_BT = K_stack{i-1}(4:6,1:3); K1_BB = K_stack{i-1}(4:6,4:6);
        
        Uval_bot_prelayer = K1_TB^(-1)*( Sval_L1top - K1_TT*Uval_L1top );  % Uval_bot_L{i-1}
        Sval_bot_prelayer = K1_BT*Uval_L1top + K1_BB*Uval_bot_prelayer;    % Sval_bot_L{i-1}
        
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