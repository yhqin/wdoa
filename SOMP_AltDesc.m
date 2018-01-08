function [Delta,X] = SOMP_AltDesc(A,B,Y,k,algType,maxIters,errTol)
%   SOMP_AltDesc - Alternating Descent algorithm usign SOMP algorithm.
%
%   This algorithm solves the nonconvex optimization problem:
%
%    min_{X,delta} ||(A + B diag(delta)) X - Y||_F + ||delta||_2 + lambda*||X||_1,2
%   
%   by a suboptimal iterative two step approach alternating between 
%   estimates of X and delta. At the first step, for a fixed delta
%   SOMP_AltDesc estimates the k nonzero entries of the jointly sparse
%   matrix X:
% 
%    min_{X} ||(A + B diag(delta)) X - Y||_F   s.t.  ||X||_0,2 <= k
%
%   Simultaneous Orthogonal Matching Pursuit (SOMP) is used to approximate
%   the solution to the above recovery problem with stopping condition the 
%   number of iterations to be as many as the sparsity level k. At the
%   second step of the algorithm delta can be estimated using the least 
%   squares inversion resulting in the SOMP_LS algorithm or the SRTLS 
%   update rule as described in [1]. The algorithm iterates between
%   estimate updates of X and delta until the norm of the current estimate 
%   and the previous iteration estimate of delta falls below the predefined
%   threshold errTol or the maximum number of iterations maxiIters has been 
%   reached.
%
%   Input arguments:
%       
%       A        - m x n dictionary (with atoms drawn from a parametric 
%                  function)
%       B        - m x n dictionary (with atoms the derivatives of those in 
%                  A)
%       Y        - m x L matrix of observations
%       k        - cardinality of the support set of X (sparsity level)
%       algType  - 'SOMP-LS' for the SOMP_LS algorithm or 'SOMP-TLS' for 
%                   the SOMP_TLS algorithm                    
%       maxIters - maximum number of iterations (optional,default is 20)
%       errTol   - error tolerance (optional,default is 1e-3)
%
%
%   Output arguments:
%
%       Delta    - vector containing the difference between the assumed 
%                  model and the estimated one (due to mismatches)
%       X        - estimate of jointly sparse solution matrix
%
%
%   Reference:
%
%       [1]A. Gretsistas & M.D. Plumbley, "An Alternating Descent algorithm
%       for the off-grid DOA estimation problem with sparsity constraints" 
%       EUSIPCO 2012.
%
%   Created: January, 2012
%   Author: A. Gretsistas
%

if nargin < 7, errTol = 1e-3; end
if nargin < 6, maxIters = 20; end
if nargin < 5, algType = 'SOMP-LS'; end
   
if nargin < 4 || isempty(A) || isempty(B) || isempty(Y)|| isempty(k) 
   error('At least four arguments are required.');
end 

[m,n] = size(A); L = size(Y,2);
A0 = A; itr = 0; Delta = zeros(n,1);
while  itr < maxIters
    itr = itr + 1;
    Delta_old = Delta;
    X = SOMP(A,Y,k);
    
    if (strcmp(algType,'SOMP-LS'))
        B_dx = zeros(m*L,n); R = zeros(m*L,1);
        for jj = 1 : L
            B_dx((jj-1)*m+1:jj*m,:) = B*diag(X(:,jj));
            R((jj-1)*m+1:jj*m,1) = Y(:,jj)-A0*X(:,jj);
        end
        Delta = real(pinv(B_dx)*R);
    elseif (strcmp(algType,'SOMP-TLS'))
        Delta = real(diag(pinv(B)*(Y-A0*X)*X'/(eye(n)+X*X')));
    else
        error('Invalid input: input argument algType wrongly inserted.');
    end
    
    A = A0 + B*diag(real(Delta));
    resnorm = norm(Delta-Delta_old);
    if (resnorm <= errTol)
        break;
    end
end
Delta = diag(Delta);