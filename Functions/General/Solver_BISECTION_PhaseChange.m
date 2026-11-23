function all_roots_x = Solver_BISECTION_PhaseChange(ComplexFun_x, x_min, x_max, dx, tol)
%The phase version of bisection method for solving a nonlinear equation ComplexFun(x)=0 within the interval [x_min, x_max],
% where ComplexFun is complex-valued while x is real-valued.
% Inputs:
% ComplexFun_x = the function handle representing ComplexFun(x),
% x_min        = the left  end point of the interval [x_min, x_max],
% x_max        = the right end point of the interval [x_min, x_max],
% dx           = the step size during searching roots,
% tol          = the tolerance value of the searched roots.
% Output:
% all_roots_x  = the roots that make ComplexFun_x(x)=0.
%
% Author: Shuanglin GUO @PIMM of ENSAM
% February 18, 2020
% Email: shuanglin.guo@ensam.eu

if nargin<4
    dx = (x_max-x_min)/10;
end
if nargin<5
    tol = 1e-5;  % the default tolerance for finding roots
end
if x_min > x_max
    error('x_min cannot be larger than x_max!');
end

delta = 5;   % unit degree

x_range = x_min:dx:x_max;
if x_range(end)<x_max
    x_range(end+1) = x_max;
end

phase_range = zeros(size(x_range));
for i = 1:length(x_range)
    x = x_range(i);
    phase_range(i) = CorrectPhase( ComplexFun_x(x) );
end

all_roots_x = [];
for i = 1:length(x_range)-1
    x1 = x_range(i);   phase1 = phase_range(i);
    x2 = x_range(i+1); phase2 = phase_range(i+1);
    
    if isnan(phase1) || isnan(phase2) || isinf(phase1) || isinf(phase2)
        root_temp = [];    % phase1 and phase2 are abnormal, there is no root in [x1, x2]
    elseif abs(phase1-phase2)>=180-delta && abs(phase1-phase2)<=180+delta          % there is a root in [x1, x2]
        flag = 0;
        while abs(x1-x2)>tol
            x0 = (x1+x2)/2; phase0 = CorrectPhase( ComplexFun_x(x0) );
            if abs(phase1-phase0)>=180-delta && abs(phase1-phase0)<=180+delta      % there is a root in [x1, x0]
                %x1 = x1; phase1 = phase1;
                x2 = x0; phase2 = phase0;
            elseif abs(phase0-phase2)>=180-delta && abs(phase0-phase2)<=180+delta  % there is a root in [x0, x2]
                x1 = x0; phase1 = phase0;
                %x2 = x2; phase2 = phase2;
            else
                flag = 1;
                break;
            end
        end
        if flag==0
            root_temp = (x1+x2)/2;
        else
            root_temp = [];
        end
    else
        root_temp = [];
    end
    
    if ~isempty(root_temp)
        all_roots_x(end+1,1) = root_temp;
    end
end

end