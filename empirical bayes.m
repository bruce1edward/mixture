rng('default')
addpath('..\functions')
iters = 20;iter = 400; % number of iteration for EM alogrithm
mcmc_steps = 8000; burn_steps = 4000; 
n = 1000; d= 20; K = floor(n*0.3); sigma = 1; b = 3;
[X,Y,Y_P,beta,Pi,inv_Pi] = generate_distribution_sparse(n, d, K, sigma, b);
X_wo = X(:,2:end);X_centered = X_wo - repmat(mean(X_wo,1), [n 1]);
%Naive 
beta_naive = X\Y_P;
%Oracle
beta_oracle = X\Y;
%Robust
beta_robust = robustfit(X_wo,Y_P,'huber');
%Mixture
Y_P1 = Y_P - mean(Y_P);control = "robust";
[beta_mixture, alpha_mixture, sigma_mixture] = fit_mixture(X_centered, Y_P1, control);
beta_mixture = [mean(Y_P - X_wo*beta_mixture) ; beta_mixture];
%EM
order = 1:n;
[beta_EM, sigma_EM] = EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, 0, beta_naive, order);
%EMM
order = 1:n;
theta = Choose_theta(n,K);
[beta_EMM, sigma_EMM] =  EM_mal_tricks(Y_P, X, iter, mcmc_steps, burn_steps, theta, beta_naive, order);


function [B_hat2, sigma_sq] = EM_mal_tricks_EB(Y_permuted, X, iters, mcmc_steps, burn_steps, theta, beta_start, order)
n = numel(order);
B_hat2 = beta_start; % initial OLS estimator of B
sigma_sq = norm(Y_permuted - X*B_hat2)^2/n; % initial OLS estimator of sigma^2
for k = 1 : iters
        Y_hat = X*B_hat2/sigma_sq; % Computing the estimator of Y_hat
        [hat_Pi,order] = mcmc_mex_mal(Y_permuted, Y_hat, order, mcmc_steps, burn_steps, theta);\
        Pi_Y = hat_Pi*Y_permuted;
        B_hat2 = X\(hat_Pi*Y_permuted);    % Get the Least square solution 
        sigma_sq = norm(Pi_Y - X*B_hat2)^2/n;
        fun = @(gamma) log(poisscdf(n-K,exp(gamma))/poisscdf(n,exp(gamma)-1)) + gamma*trace(hat_Pi, eye(n))
        theta = fminbnd(fun,x1,x2);
end
end