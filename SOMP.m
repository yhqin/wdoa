function [X,actSet] = SOMP(A,Y,maxIters,errTol)
%   SOMP - Simultaneous Orthogonal Matching Pursuit.
%
%   Input arguments:
%       
%       A        - m x n dictionary
%       Y        - m x L matrix of observations
%       maxIters - maximum number of iterations (optional,default is m)
%       errTol   - error tolerance (optional,default is 1e-5)
%
%
%   Output arguments:
%
%       X        - estimate of jointly sparse solution matrix
%       actSet   - array of support set indices
%
%
%   Created: January, 2012
%   Author: A. Gretsistas
%

if nargin < 4, errTol = 1e-5; end
if nargin < 3, maxIters = size(A,1); end

if nargin < 2 || isempty(A) || isempty(Y)
   error('At least two arguments are required.');
end 

L = size(Y,2);
n = size(A,2);
lambdareg = 1e-8;

R = Y;

actSet = [];
X = zeros(n,L);
itr = 1;

while (itr <= maxIters) && (norm(R,'fro') > errTol)
    [~,indx] = max(sum(abs(A'*R),2));
    actSet = [actSet indx];
    Phi = A(:,actSet);
    regI = lambdareg*eye(itr);
    for j = 1 : L
        X(actSet,j) = (Phi'*Phi + regI)\(Phi'*Y(:,j));
    end
    Yhat = Phi*X(actSet,:);
    R = Y - Yhat;
    itr = itr + 1;
end