function [ROOT_x, num_fun_eval] = Solver_MRCM1D(Fun_x, x_min, x_max, dx)
%MRCM1D Find zeros of a univariate function by modulus-ratio convergence.
%
% ROOT_x = MRCM1D(Fun_x, x_min, x_max, dx)
% [ROOT_x, num_fun_eval] = MRCM1D(...)
%
% The second output is the number of scalar calls to Fun_x. Function values
% at grid points shared by adjacent scanning intervals and at the inherited
% endpoints/centre of an interval-refinement step are reused instead of being
% evaluated again.
%
% Inputs:
% Fun_x = function handle to represent Fun(x).
% x_min = the left  end point in the range [x_min, x_max].
% x_max = the right end point in the range [x_min, x_max].
% dx = the interval step.
% Outputs:
% ROOT_x = the roots of the equation Fun(x)=0 found within the range [x_min, x_max].
%
% Author: Shuanglin GUO @PIMM of ENSAM
% February 19, 2020
% shuanglin.guo@ensam.eu
% 
% Reference:
% Feng Zhu, Bin Wang & Zheng-hua Qian. A numerical algorithm to solve
% multivariate transcendental equation sets in complex domain and its
% application in wave dispersion curve characterization. Acta Mechanica,
% volume 230, pages 1303-1321 (2019).

step_node = dx/5;

num_fun_eval = 0;
ROOT_x = [];

% Algorithm 2.1.1 of the reference: seven nodes per scanning interval, of which
% the last two are also the first two nodes of the next scanning interval.
num_main_interval = numel(x_min:dx:x_max);
if num_main_interval == 0
    return;
end

num_main_node = 7;
ROOT_x = zeros(1, num_main_interval);
num_root = 0;
x_value = zeros(1, num_main_node);
abs_value = zeros(1, num_main_node);

for k = 1:num_main_interval
    if k == 1
        x_value = x_min + (0:num_main_node-1)*step_node;
        eval_idx = 1:num_main_node;
    else
        % Inherit the two overlapping nodes from the previous interval.
        x_value(1:2) = x_value(6:7);
        abs_value(1:2) = abs_value(6:7);

        % Evaluate only the five new nodes.
        x_value(3:7) = x_value(2) + (1:5)*step_node;
        eval_idx = 3:7;
    end

    for j = eval_idx
        abs_value(j) = abs(Fun_x(x_value(j)));
        num_fun_eval = num_fun_eval + 1;
    end

    absFun_min = Inf;
    idx_candidate = [];
    for j = 1:num_main_node
        if abs_value(j) < absFun_min
            absFun_min = abs_value(j);
            idx_candidate = j;
        end
    end

    % A local minimum must be strictly inside the scanning interval.
    if ~isempty(idx_candidate) && idx_candidate > 1 && idx_candidate < num_main_node
        [root_x, num_local_eval] = JudgeNullPoints_1D(...
            Fun_x, x_value(idx_candidate), step_node,...
            abs_value(idx_candidate-1),...
            abs_value(idx_candidate),...
            abs_value(idx_candidate+1));
        num_fun_eval = num_fun_eval + num_local_eval;

        if ~isempty(root_x)
            num_root = num_root + 1;
            ROOT_x(num_root) = root_x;
        end
    end
end

if num_root == 0
    ROOT_x = [];
else
    ROOT_x = ROOT_x(1:num_root);
    ROOT_x(ROOT_x > x_max) = [];
end
end


function [root_x, num_fun_eval] = JudgeNullPoints_1D(...
    Fun_x, x_candidate, step_node, absFun_left, absFun_center, absFun_right)
% Algorithm 2.1.2 of the reference with inherited endpoint and centre values.
% In every refinement step the new interval is centred at x_candidate and
% its endpoints are x_candidate +/- step_node. These three points were already
% evaluated in the preceding step, so only the other 18 nodes are evaluated.

M = 1e3;
num_iter = 16;
num_node = 21;

num_fun_eval = 0;
absFun_min = Inf;
x_min_candidate = [];
idx_candidate = [];
absFun_min_record = M*ones(num_iter, 1);

% First refined interval. Its endpoints and centre are known from the main scanning interval.
node_step = step_node/10;
x_value = x_candidate + (-10:10)*node_step;
abs_value = nan(1, num_node);
known_idx = [1, 11, 21];
abs_value(known_idx) = [absFun_left, absFun_center, absFun_right];
eval_idx = [2:10, 12:20];

for k = 1:num_iter
    if k > 1
        % The endpoints and centre of the new interval are the three
        % consecutive nodes around the previous minimum.
        previous_abs_value = abs_value;
        previous_idx_candidate = idx_candidate;

        node_step = node_step/10;
        x_value = x_min_candidate + (-10:10)*node_step;
        abs_value = nan(1, num_node);
        abs_value([1, 11, 21]) = previous_abs_value(...
            previous_idx_candidate-1:previous_idx_candidate+1);
        eval_idx = [2:10, 12:20];
    end

    for j = eval_idx
        abs_value(j) = abs(Fun_x(x_value(j)));
        num_fun_eval = num_fun_eval + 1;
    end

    % Keep the original strict comparison order and tie handling.
    for j = 1:num_node
        if abs_value(j) < absFun_min
            absFun_min = abs_value(j);
            x_min_candidate = x_value(j);
            idx_candidate = j;
        end
    end

    absFun_min_record(k) = absFun_min;
end

if absFun_min_record(end)*M < absFun_min_record(1)
    root_x = x_min_candidate;
else
    root_x = [];
end
end