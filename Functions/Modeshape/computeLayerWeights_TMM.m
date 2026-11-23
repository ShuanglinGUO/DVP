function weights_layer = computeLayerWeights_TMM(k3_layer, U_layer, D_layer, d_layer, fq, cp, method_LinearSolver)
%Compute layer weights of partial waves using TMM method.
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

num_layer = length(k3_layer);
xi = 2*pi*fq/cp;

t_layer = cell(1,num_layer);
T_stack = cell(1,num_layer);
for i = 1:num_layer
    k3 = k3_layer{i}; U = U_layer{i}; D = D_layer{i}; d = d_layer(i);
    
    T1 = [U; D];
    T2 = diag( exp(1i*k3*d) );
    t_layer{i} = T1*T2*T1^(-1);
    
    % the recursive process
    if i==1
        T_stack{1} = t_layer{1};
    else
        T_stack{i} = t_layer{i}*T_stack{i-1};
    end
end

Uval_L1top = solveLinearHomoEqSet(T_stack{end}(4:6,1:3), method_LinearSolver);
Sval_L1top = [0; 0; 0];    % sigma33(x3=0)=0, sigma23(x3=0)=0, sigma13(x3=0)=0
USval_L1top = [Uval_L1top; Sval_L1top];

USval_top_L = cell(1,num_layer);
for i = 1:num_layer
    if i==1
        USval_top_L{1} = USval_L1top;
    else
        USval_bot_prelayer = t_layer{i-1}*USval_top_L{i-1};  % USval_bot_L{i-1}
        USval_top_L{i} = USval_bot_prelayer;  % continuity condition: USval_top_L{i} = USval_bot_L{i-1}
    end
end

weights_layer = cell(1,num_layer);
for i = 1:num_layer
    U = U_layer{i};
    D = D_layer{i};
    weights_layer{i} = [U; D]^(-1)*USval_top_L{i};
end

end