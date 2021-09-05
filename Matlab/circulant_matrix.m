% This function generates a lower triangular Toeplitz square matrix of
% dimension (n + nu) whose first column is equal to the (nu + 1) impulse
% response co-efficients of channel H. n is the dimension of each sub-block
% i.e. each sub-block has n symbols and nu is the memory of channel
% H is a row vector of length nu

function X = lower_toeplitz(n, nu, H)

% n = 8;
% nu = 2;
% H = [1 2 3];
X(:,1) = [H zeros(1, n+nu-length(H))].';
    
for i = 2 : n+nu
    X(:,i) = circshift(X(:,i-1),1);
end
