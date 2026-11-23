function [k3, U, D] = solveStrohMatrix(C, rho, fq, cp)
%Stroh formalism to compute the eigen value k3 and eigen vectors U and D.
% k3 is the wavenumber in x3 direction.
% U and D are the partial wave amplitude matrix of displacement and stress.
%
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

factor_scale = 1;    % scale factor

omega = 2*pi*fq;
k1 = omega/cp;

D11 = [C(1,1), C(1,6), C(1,5); C(1,6), C(6,6), C(5,6); C(1,5), C(5,6), C(5,5)];
D33 = [C(5,5), C(4,5), C(3,5); C(4,5), C(4,4), C(3,4); C(3,5), C(3,4), C(3,3)];
D13 = [C(1,5), C(1,4), C(1,3); C(5,6), C(4,6), C(3,6); C(5,5), C(4,5), C(3,5)];
D31 = [C(1,5), C(5,6), C(5,5); C(1,4), C(4,6), C(4,5); C(1,3), C(3,6), C(3,5)];

% Stroh matrix of the coupling waves, the eigen value is k3
S = [-k1*D33^(-1)*D31, -1i*D33^(-1);
    -1i*( k1^2*(D11-D13*D33^(-1)*D31) - rho*omega^2*eye(3) ), -k1*D13*D33^(-1)];
[UD, k3] = eig( S ); k3 = diag(k3);

[~, idx] = sort(real(k3)); k3 = k3(idx); UD = UD(:,idx);

U = UD(1:3,:);               % amplitude of displacement u1, u2, u3
D = UD(4:6,:)/factor_scale;  % amplitude of stress sigma31, sigma32, sigma33

end