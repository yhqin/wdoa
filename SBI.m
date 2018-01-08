function res = SBI(paras)

% res = SBI(paras)
% 
% SBI(paras) performs DOA estimation using Sparse Bayesian Inference
% 
% Input:
% paras.Y: M * T matrix, sensor measurements at all snapshots
% paras.A: M * N matrix, columns are the steering vectors for different directions
% paras.B: M * N matrix, columns are derivatives of the steering vectors wrt. different directions
% paras.sigma2: initialization of noise variance
% paras.alpha: initialization of alpha
% paras.beta: initialization of beta
% paras.rho: rho
% paras.resolution: grid resolution for the directions
% paras.maxiter: maximum iteration
% paras.tol: stopping criterion
% paras.isKnownNoiseVar: true if known variance, false if unknown
% paras.K: number of sources
% 
% Output:
% res.mu: mean estimation
% res.Sigma: variance estimation
% res.sigma2: estimated noise variance
% res.sigma2seq: estimated noise variance at all iterations
% res.alpha: reconstructed alpha
% res.beta: reconstructed beta
% res.iter: iteration used in the algorithm
% res.ML: maximum likelihood function value at all iterations
% 
% Written by Zai Yang, 19 Jul, 2011
% reference: Zai Yang, Lihua Xie, and Cishen Zhang, 
%    "Off-grid direction of arrival estimation using sparse Bayesian inference"

eps = 1e-16;

Y = paras.Y;
A = paras.A;
B = paras.B;

[M, T] = size(Y);
N = size(A, 2);

alpha0 = 1 / paras.sigma2;
rho = paras.rho / T;
beta = paras.beta;
alpha = paras.alpha;
r = paras.resolution;

maxiter = paras.maxiter;
tol = paras.tolerance;

if isfield(paras, 'isKnownNoiseVar') && ~isempty(paras.isKnownNoiseVar)
    isKnownNoiseVar = paras.isKnownNoiseVar;
else
    isKnownNoiseVar = false;
end

if isKnownNoiseVar
    a = 1;
    b = T * M * paras.knownsigma2;
else
    a = 1e-4;
    b = 1e-4;
end

if isfield(paras, 'K') && ~isempty(paras.K)
    K = paras.K;
else
    K = min(T, M-1);
end

idx = [];
BHB = B' * B;
converged = false;
iter_beta = 1;
iter = 0;
ML = zeros(maxiter,1);
alpha0seq = zeros(maxiter,1);

while ~converged
    iter = iter + 1;
    
    Phi = A;
    Phi(:,idx) = A(:,idx) + B(:,idx) * diag(beta(idx));
    
    alpha_last = alpha;
    
    C = 1 / alpha0 * eye(M) + Phi * diag(alpha) * Phi';
%     Sigma = diag(alpha) - diag(alpha) * Phi' / C * Phi * diag(alpha);
    Cinv = inv(C);
    Sigma = diag(alpha) - diag(alpha) * Phi' * Cinv * Phi * diag(alpha);
    mu = alpha0 * Sigma * Phi' * Y;
    
    
    gamma1 = 1 - real(diag(Sigma)) ./ (alpha + eps);
    
    % update alpha
    musq = mean(abs(mu).^2, 2);
    
    alpha = musq + real(diag(Sigma));
    if rho ~= 0
        alpha = -.5 / rho + sqrt(.25 / rho^2 + alpha / rho);
    end
    
    
    % update alpha0
    resid = Y - Phi * mu;
    alpha0 = (T * M + a - 1) / (norm(resid, 'fro')^2 + T / alpha0 * sum(gamma1) + b);
    alpha0seq(iter) = alpha0;
    
    % stopping criteria
    if norm(alpha - alpha_last)/norm(alpha_last) < tol || iter >= maxiter
        converged = true;
        iter_beta = 5;
    end
    
    temp = 0;
    for t = 1:T
        temp = temp + real(Y(:,t)' * Cinv * Y(:,t));
    end
    ML(iter) = -T * real(log(det(C))) - temp + (a-1) * log(alpha0) - b * alpha0 - rho * sum(alpha);

    
    % update beta
    [temp, idx] = sort(alpha, 'descend');
    idx = idx(1:K);

%     [peaks, idx] = findpeaks(alpha,'sortstr','descend');
%     if length(idx) > K
%         idx = idx(1:K);
%     end
    temp = beta;
    beta = zeros(N,1);
    beta(idx) = temp(idx);
      
    P = real(conj(BHB(idx,idx)) .* (mu(idx,:) * mu(idx,:)' + T * Sigma(idx,idx)));
    v = zeros(length(idx), 1);
    for t = 1:T
        v = v + real(conj(mu(idx,t)) .* (B(:,idx)' * (Y(:,t) - A * mu(:,t))));
    end
    v = v - T * real(diag(B(:,idx)' * A * Sigma(:,idx)));
    temp1 = P \ v;
    if any(abs(temp1) > r/2) || any(diag(P) == 0)
        for i = 1:iter_beta
            for n = 1:K
                temp_beta = beta(idx);
                temp_beta(n) = 0;
                beta(idx(n)) = (v(n) - P(n,:) * temp_beta) / P(n,n);
                if beta(idx(n)) > r/2
                    beta(idx(n)) = r/2;
                end
                if beta(idx(n)) < -r/2
                    beta(idx(n)) = -r/2;
                end
                if P(n,n) == 0
                    beta(idx(n)) = 0;
                end
            end
        end
    else
        beta = zeros(N,1);
        beta(idx) = temp1;
    end  
    
end

res.mu = mu;
res.Sigma = Sigma;
res.beta = beta;
res.alpha = alpha;
res.iter = iter;
res.ML = ML(1:iter);
res.sigma2 = 1/alpha0;
res.sigma2seq = 1./alpha0seq(1:iter);

end
