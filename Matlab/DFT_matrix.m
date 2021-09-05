
% This function generates the orthonormal Discrete Fourier Transform (DFT)
% matrix Q whose (l,k) element is given by 
% Q(l,k) = 1/sqrt(n) * exp (-j * 2 * pi * l * k) where 0<=l,k<=n-1
% n = block length
% dimension of Q is n by n

function Q = DFT_matrix(n)

for l=1:n
    for k=1:n
        Q(l,k) = sqrt(1/n) * exp (-j * (2 * pi * (l-1) * (k-1))/n);
    end
end