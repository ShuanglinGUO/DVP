function C = StiffnessMatrixRotation(C_star, phi)
%Compute stiffness matrix rotation according to the tensor formula.
% Note: phi = the rotation angle which should be with the unit: rad
% Author: Shuanglin GUO @PIMM of ENSAM
% June 10, 2018
% Email: shuanglin.guo@ensam.eu

% the rotation matrix
beta = [cos(phi), sin(phi), 0;
       -sin(phi), cos(phi), 0;
               0,        0, 1];

C = zeros(6,6);
for m = 1:3
    for n = 1:3
        Idx1 = TensorIndex2MatrixIndex(m,n);
        for o = 1:3
            for p = 1:3
                Idx2 = TensorIndex2MatrixIndex(o,p);
                
                sum_temp = 0;
                for i = 1:3
                    for j = 1:3
                        idx1 = TensorIndex2MatrixIndex(i,j);
                        for k = 1:3
                            for l = 1:3
                                idx2 = TensorIndex2MatrixIndex(k,l);
                                
                                sum_temp = sum_temp + beta(m,i)*beta(n,j)*beta(o,k)*beta(p,l)*C_star(idx1,idx2);
                            end
                        end
                    end
                end
                
                C(Idx1,Idx2) = sum_temp;
            end
        end
    end
end
end


function idx = TensorIndex2MatrixIndex(i,j)

if i==1 && j==1
    idx = 1;
elseif i==2 && j==2
    idx = 2;
elseif i==3 && j==3
    idx = 3;
elseif (i==2 && j==3) || (i==3 && j==2)
    idx = 4;
elseif (i==1 && j==3) || (i==3 && j==1)
    idx = 5;
elseif (i==1 && j==2) || (i==2 && j==1)
    idx = 6;
else
    error('There is something wrong!');
end

end