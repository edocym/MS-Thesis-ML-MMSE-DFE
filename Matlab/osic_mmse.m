% This function uses the osic mmse algorithm to estimate the transmitted
% symbol
% H_eqv = the channel
% Y_eqv = received column vector
% s = transmitted symbols
% snr = SNR
% block_length = # of symbols transmitted in a block (not sub-block). This is
%                  needed to normalize the energy of the symbols in Z_mmse_osic 
%                  in the same way as in transmitted symbols
% Z_mmse_osic = estimated received symbols

% VERY VERY IMPORTANT each element of Z_mmse_osic should have the same
% energy as transmitted symbols

function Z_mmse_osic = osic_mmse(H_eqv, Y_eqv, snr, block_length)

[Q,R] = qr(H_eqv);
D = diag(R).^2;
        
for kk = 1:length(Y_eqv)
    G = inv(H_eqv' * H_eqv + (1/snr) * eye(size(H_eqv,2))) * H_eqv';
%     G = inv(H_eqv' * H_eqv + (1/snr) * eye(size(H_eqv,2)));
    [val,ind] = max(D);
    Y_ind = G(ind,:) * Y_eqv;
    received_bits = [sign(real(Y_ind))+1 sign(imag(Y_ind))+1]./2;
    Z_mmse_osic(ind,:) = (qpsk_from_bits(received_bits) / sqrt(block_length));
    Y_eqv = Y_eqv - H_eqv(:,ind) * Z_mmse_osic(ind,:);
    H_eqv = H_eqv - [zeros(size(H_eqv,1),ind-1) ...
                             H_eqv(:,ind) zeros(size(H_eqv,1),size(H_eqv,2)-ind)];
    D(ind) = 0;
end