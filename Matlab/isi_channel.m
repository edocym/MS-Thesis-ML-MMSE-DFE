%channel = [randn(m,n) + j * randn(m,n)] returns a m * n matrix of random 
%numbers chosen from the normal distribution. This function is creating a
%channel with memory nu i.e. a nu tap channel whose co-efficients
%are normally distributed random numbers.

function h = isi_channel(nu);

sigma_sq = rand(1,nu+1);
sigma_sq = sigma_sq ./sum(sigma_sq);

h = sqrt(sigma_sq) .* [randn(1,nu+1) + j * randn(1, nu+1)];