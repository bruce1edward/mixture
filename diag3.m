% Sparse Permutation only for k = 4
n = 10;
d = 1;
k = 5;
noise = 0.1;
pi = perms(1:n);
order = 1:n;
spi = pi(sum(pi ~= order,2) <= k,:);
[nsp,d1] = size(spi);
Pi = eye(n);

%Exact Expectation
SP = spi(unidrnd(nsp),:);
[x,y,coefs] = generate_distribution(n, d, noise, SP);
SNR = coefs^2/noise;
p_density = zeros(nsp,1);
hat_SPi = zeros(n,n);
for j = 1 : nsp 
p_density(j) = dot(y,x(spi(j,:),:)*coefs)/noise;
end
p_density =  p_density - max(p_density);
p_density = softmax(p_density);

for l = 1 : nsp
 hat_SPi = hat_SPi + Pi(spi(l,:),:)*p_density(l);
end

mcmc_steps = 4*10^6;
burn_steps = floor(mcmc_steps/4); % 1/4 burn in 
% Linear indexing
index = 1:n;
index1 = @(order) (index - 1)*n + order; % transfering double index to linaer index

%Local Sampling
order_tmm = zeros(mcmc_steps,1);
S = generate_randerang(n, k);
[S_C,order] = generate_off_Pi(S,n);
S_abs = length(S);
p_prior = [0.4 0.4]; % P(Keep) = 0.4/ P(Add) = 0.4 / P(Delete) = 0.2 
for m_step = 1: mcmc_steps
    p1 = dot(y, x(order,:)*coefs)/noise;
    %Sample 
    u = rand;
    if S_abs == 5
       %keep or delete
       if u <= p_prior(1)
        %keep
        [S2,S_C2,order2,S_abs2] = keep(S,S_C,order,S_abs);
       end
       if p_prior(1) < u <= 1
        %delete
        [S2,S_C2,order2,S_abs2] = del(S,S_C,order,S_abs,n);
       end
    end
   
    if S_abs == 3 || S_abs == 2 || S_abs == 4
        %delete, add or keep 
       if u > p_prior(1) + p_prior(2) 
        %delete
        [S2,S_C2,order2,S_abs2] = del(S,S_C,order,S_abs,n);
       end
       if p_prior(1) < u <= p_prior(1) + p_prior(2)
        %add
        [S2,S_C2,order2,S_abs2] = add(S,S_C,order,S_abs,n);
       end
       if u <= p_prior(1) 
        %keep
        [S2,S_C2,order2,S_abs2] = keep(S,S_C,order,S_abs);
       end
    end
    
    if S_abs == 0
       %add 
       [S2,S_C2,order2,S_abs2] = add(S,S_C,order,S_abs,n);
    end
    p2 = dot(y, x(order2,:)*coefs)/noise;
    %[Lia2, Locb2] = ismember(order_, spi, 'rows');
    %if Lia2 == 1
       %order_tmm(m_step,3) = Locb2;
    %end
     if p1 < p2 
        %order_tmm(m_step,1) = 1;
        %Update the support of permutation
        order = order2;
        S = S2;
        S_C = S_C2;
        S_abs = S_abs2;
     else
        if rand < exp(p2 - p1)  %since they're log probabilities
           %order_tmm(m_step,1) = 1;
           %Update the support of permutation
           order = order2;
           S = S2;
           S_C = S_C2;
           S_abs = S_abs2;
        end
     end
     %S_t(:,m_step) = S;
     %S_C_t(:,m_step) = S_C;
    if m_step >= burn_steps
     [Lia, Locb] = ismember(order, spi, 'rows');
     if Lia == 1
        order_tmm(m_step) = Locb;
     end
    end
end
p_density_e1 = tabulate(order_tmm(10^6:4*10^6));
p_density_e1 = p_density_e1(:,3)/100;
error_p1 = p_density - p_density_e1;
norm(error_p1)

index_spi = 1:nsp;
figure
hold on 
p1 = plot(index_spi, p_density, 'b', 'LineWidth', 1.2);
xlabel('Index of sparse permutation')
ylabel('Probability')
title(['Posterior of \Pi (n = 10 and K = 4) SNR = ' num2str(SNR)])
hold off

figure
hold on 
p1 = plot(index_spi, p_density_e1, 'b', 'LineWidth', 1.2);
xlabel('Index of sparse permutation')
ylabel('Probability')
title(['Posterior of \Pi (n = 10 and K = 4) SNR = ' num2str(SNR)])
hold off

figure
plot(error_p1) %  
xlabel('Index of sparse permutation')
ylabel('True - Empirical probability based on local sampling for 4 - sparse')
