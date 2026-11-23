function vec = solveLinearHomoEqSet(A, method, index_FreeVariable)
%Solve a linear homogeneous equation set, [A]{vec}={0}, where A should be a square matrix with dimension n x n.
% Author: Shuanglin GUO @PIMM of ENSAM
% April 12, 2024
% Email: shuanglin.guo@ensam.eu

if nargin<2
    method = 'indeterminacy';
end
if size(A,1) ~= size(A,2)
    error('Matrix A must be a square matrix.');
else
    n = size(A,1);
end

switch method
    case 'indeterminacy'  % one degree on indeterminacy
        % take the first variable as the free variable, this method is based on the least square approach
        vec = [1; -A(:,2:end)\A(:,1)];
        
    case 'indeterminacy_icol'  % one degree on indeterminacy
        % take the i-th variable (transferred from exterior) as the free variable
        if ~exist('index_FreeVariable', 'var')
            error('No index_FreeVariable input!');
        else
            i = index_FreeVariable;
        end
        vec_right = A(:,i);       % take the i-th column of matrix A
        matrix_left = A(:,[1:i-1,i+1:end]); vec_temp = -matrix_left\vec_right;
        vec = zeros(n,1);
        vec(i) = 1;
        vec([ 1:i-1, i+1:end ]) = vec_temp;
        
    case 'normalization'
        sub_A = A(1:n-1,1:n-1); vec_n = A(1:n-1,n);
        
        Delta = zeros(n,1);
        for i = 1:n
            if i==n
                Delta(n) = det( sub_A );
            else
                sub_Ai = sub_A;
                sub_Ai(:,i) = vec_n;
                Delta(i) = -det( sub_Ai );
            end
        end
        
        X = 1/norm( Delta );
        vec = Delta*X;
        
    case 'null'  % use MATLAB built-in function, null
        if rank(A)==n-1
            vec = null(A);
        else
            error('A is not a matrix having one degree on indeterminacy, please use other method!');
        end
        
    case 'eig'  % eigenvalue decomposition-based method
        if issparse(A)
            [eig_vector, eig_value] = eigs(A);   % if A is defined as a sparse matrix, eigs should be used
        else
            [eig_vector, eig_value] = eig(A);
        end
        [~, idx] = min(abs( diag(eig_value) ));  % as A is singular, its minimum absolute eigen value should approach to zero
        vec = eig_vector(:,idx);                 % the eigen vector corresponding to this eigen value should be the vector wanted
        
    case 'svd'  % singular value decomposition-based method
        [~, S, V] = svd(A);
        [~, idx_min] = min(abs( diag(S) ));
        vec = V(:,idx_min);
        
    otherwise
        error('There are only six cases for method: ''indeterminacy'', ''indeterminacy_icol'', ''normalization'', ''null'', ''eig'' and ''svd''.');
end

end