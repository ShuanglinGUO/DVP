function [fq, cg] = computeGroupVelocity_Isotropic(fq, cp, options)
%Compute the dispersion curve of group velocity in isotropic material based on the normal formula:
% cg = partial derivative of omega with respect to xi (do not consider the anisotropic effect).
% Note #1: the current function treats a single mode at a time.
% Note #2: the curve (fq, cp) should be sorted externally.
%
% Author: Shuanglin GUO @PIMM of ENSAM
% July 4, 2020
% Email: shuanglin.guo@ensam.eu

% set default parameters for options---------------------------------------
options_default.tol = 1e-5;
options_default.num_sampling = 5000;
options_default.downsampling = 'OFF';  % 2 options: 'ON' and 'OFF'
options_default.sampling4what = 'FQ';  % 2 options: 'FQ' and 'CP'
if nargin<3
    options = options_default;
else
    options = setOptions(options, options_default);
end
tol = options.tol;
num_sampling = options.num_sampling;
downsampling = options.downsampling;
sampling4what = options.sampling4what;
% end to set default parameters for options--------------------------------


if length(fq)~=length(cp)
    error('The length of fq should equal to the length of cp!');
end

% delete the repeated points
for i = length(fq):-1:2
    if abs( (fq(i)-fq(i-1))/(fq(i)+eps) ) + abs( (cp(i)-cp(i-1))/(cp(i)+eps) ) < tol
        fq(i) = []; cp(i) = [];
    end
end

% up-sampling to make the computed cg curve more smooth
switch sampling4what
    case 'FQ'  % sampling for frequency (FQ), if cp has monotonicity with respect to fq, use this case
        fq_sp = linspace( min(fq), max(fq), num_sampling );
        cp_sp = spline( fq, cp, fq_sp );
    case 'CP'  % sampling for phase velocity (CP), if there exists negative group velocity, use this case
        cp_sp = linspace( min(cp), max(cp), num_sampling );
        fq_sp = spline( cp, fq, cp_sp );
    otherwise
        error('sampling4what only has two options: ''FQ'' and ''CP''!');
end

cg_sp = zeros(size(cp_sp));
for i = 1:length(fq_sp)-1
    Dcp_Dfq = ( cp_sp(i)-cp_sp(i+1) )/( fq_sp(i)-fq_sp(i+1) );
    cg_sp(i) = cp_sp(i)^2/( cp_sp(i) - fq_sp(i)*Dcp_Dfq );
end
cg_sp(end) = cg_sp(end-1);

% down-sampling to match the input frequency, fq, if downsampling='ON'
if strcmp(downsampling, 'ON')
    cg = spline( fq_sp, cg_sp, fq );
elseif strcmp(downsampling, 'OFF')
    cg = cg_sp; fq = fq_sp;
else
    error('downsampling has only two options: ''ON'' and ''OFF''!');
end

end