function [y, t_layer, T_stack] = DispersionEQ_TMM(LayerPro, fq, cp)
%Construct the dispersion equation using transfer matrix method (TMM).
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

num_layer = length(LayerPro);
% k1 = 2*pi*fq/cp;  % wavenumber in x1 direction

t_layer = cell(1,num_layer);
T_stack = cell(1,num_layer);
for i = 1:num_layer
    d = LayerPro(i).d;
    [k3, U, D] = solveStrohMatrix(LayerPro(i).rotC, LayerPro(i).rho, fq, cp);
    
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

y = det( T_stack{end}(4:6,1:3) );

end