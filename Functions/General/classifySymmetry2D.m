function [X_sym, Y_sym, modeshape_sym,...
    X_anti, Y_anti, modeshape_anti,...
    X_none, Y_none, modeshape_none] = classifySymmetry2D(...
    X_cdd, Y_cdd, modeshape_disp_cdd, options)
%Classify the candidate points (X_cdd, Y_cdd) into symmetric and anti-symmetric modes.
%If a candidate point has no any symmetry, it will be assigned into none symmetric mode.
% Author: Shuanglin GUO @PIMM of ENSAM
% June 12, 2020
% Email: shuanglin.guo@ensam.eu

% set default parameters for options
options_default.tol = 1e-3;
if nargin<4
    options = options_default;
else
    options = setOptions(options, options_default);
end
tol = options.tol;


[num_disp, num_x3] = size( modeshape_disp_cdd(:,:,1) );
if num_disp ~= 3
    error('Displacement must contain u1, u2 and u3 components in this version!');
end
if mod(num_x3,2)==0      % num_x3 is an even number
    halfnum_x3 = num_x3/2;
elseif mod(num_x3,2)==1  % num_x3 is an odd number
    halfnum_x3 = (num_x3+1)/2;
else
    error('num_x3 should be an integer!');
end

count_sym  = 0;  X_sym = [];  Y_sym = [];  modeshape_sym = [];
count_anti = 0; X_anti = []; Y_anti = []; modeshape_anti = [];
count_none = 0; X_none = []; Y_none = []; modeshape_none = [];

for i = 1:length(X_cdd)
    modeshape_u1 = modeshape_disp_cdd(1,:,i); u1_half1 = modeshape_u1(1:halfnum_x3); u1_half2 = modeshape_u1(end:-1:end-halfnum_x3+1);
    modeshape_u2 = modeshape_disp_cdd(2,:,i); u2_half1 = modeshape_u2(1:halfnum_x3); u2_half2 = modeshape_u2(end:-1:end-halfnum_x3+1);
    modeshape_u3 = modeshape_disp_cdd(3,:,i); u3_half1 = modeshape_u3(1:halfnum_x3); u3_half2 = modeshape_u3(end:-1:end-halfnum_x3+1);

    if sum(abs(u1_half1+u1_half2))<tol && sum(abs(u2_half1+u2_half2))<tol && sum(abs(u3_half1-u3_half2))<tol  % anti-symmetric mode
        count_anti = count_anti + 1;
        X_anti(count_anti) = X_cdd(i);
        Y_anti(count_anti) = Y_cdd(i);
        modeshape_anti(:,:,count_anti) = modeshape_disp_cdd(:,:,i);

    elseif sum(abs(u1_half1-u1_half2))<tol && sum(abs(u2_half1-u2_half2))<tol && sum(abs(u3_half1+u3_half2))<tol  % symmetric mode
        count_sym = count_sym + 1;
        X_sym(count_sym) = X_cdd(i);
        Y_sym(count_sym) = Y_cdd(i);
        modeshape_sym(:,:,count_sym) = modeshape_disp_cdd(:,:,i);

    else  % none symmetric mode
        count_none = count_none + 1;
        X_none(count_none) = X_cdd(i);
        Y_none(count_none) = Y_cdd(i);
        modeshape_none(:,:,count_none) = modeshape_disp_cdd(:,:,i);
    end
end

end