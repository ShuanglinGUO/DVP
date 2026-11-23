function [C_class, C_class_number] = distinguishStiffnessMatrixClass(C)
% Author: Shuanglin GUO @NDT Lab of NCHU
% May 18, 2026
% Email: guo.shuanglin@nchu.edu.cn

tol = 1;

if sum(abs([C(1,4),C(1,5),C(1,6),C(2,4),C(2,5),C(2,6),C(3,4),C(3,5),C(3,6),C(4,5),C(4,6),C(5,6)]))<tol
    if sum(abs([C(1,3)-C(1,2), C(2,3)-C(1,2), C(2,2)-C(1,1), C(3,3)-C(1,1), C(4,4)-(C(1,1)-C(1,2))/2, C(5,5)-(C(1,1)-C(1,2))/2, C(6,6)-(C(1,1)-C(1,2))/2]))<tol
        C_class = 'Cisot';  % isotropic class
        C_class_number = 6;
    elseif sum(abs([C(1,3)-C(1,2), C(2,3)-C(1,2), C(2,2)-C(1,1), C(3,3)-C(1,1), C(5,5)-C(4,4), C(6,6)-C(4,4)]))<tol
        C_class = 'Ccubi';  % cubic class
        C_class_number = 5;
    elseif sum(abs([C(1,3)-C(1,2), C(3,3)-C(2,2), C(4,4)-(C(2,2)-C(2,3))/2, C(6,6)-C(5,5)]))<tol
        C_class = 'Ctran';  % transversely isotropic class
        C_class_number = 4;
    else
        C_class = 'Corth';  % orthotropic class
        C_class_number = 3;
    end

elseif sum(abs([C(1,4),C(1,5),C(2,4),C(2,5),C(3,4),C(3,5),C(4,6),C(5,6)])) < tol...
        && sum(abs([C(1,6),C(2,6),C(3,6),C(4,5)])) > tol
    C_class = 'Cmono';      % monoclinic class
    C_class_number = 2;

else
    C_class = 'Ctric';      % triclinic class
    C_class_number = 1;
end

end