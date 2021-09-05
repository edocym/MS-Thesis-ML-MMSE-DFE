close all; clear all; clc;

% input bit energy
Ex = 1;

%input SNR in dB
SNR = [5:5:25];
snr_linear = 10.^(SNR/10);

% sub-block length is n i.e. one sub-block contains n symbols (symbols can
% be real as PAM or complex as QPSK
n = 16;
bits_per_symbol = 2;

%channel memory nu so channel will have nu+1 taps
nu = 2;

%number of transmit antenna
Mt = 2;

%changes the state of the random generator
randn('state', sum(100*clock));
rand('state', sum(100*clock));

%run = number of transmission
run = 10000;

Q = DFT_matrix(n);

for k = 1:length(snr_linear)
    snr = snr_linear(k)
    
    total_error_zfe = 0;
    total_error_mmse_le = 0;
    total_error_zf_dfe = 0;
    total_error_mmse_dfe = 0;
    total_error_mmse_osic = 0;
    
    
    for l = 1:run
    
        %generate sub-blocks of length n symbols. Ckl means l-th sub-block of tx
        %antenna k
%         C11_row = sub_block_gen(n,Mt);
%         C21_row = sub_block_gen(n,Mt);
%         C11 = transpose(C11_row);
%         C21 = transpose(C21_row);
        
        % For 2 sub-blocks there are 2n symbols
        trans_bits = randint(1,bits_per_symbol * 2 * n);
        sym = qpsk_from_bits(trans_bits);
        symbols = sym/sqrt(2 * n);
        C11 = symbols(1:n).';
        C21 = symbols(n+1:2*n).';

        %ISI CHANNEL    

        % hkl specifies the channel from Tx antenna l to Rx antenna k
        % all are IID gaussian 3 tap channel so a vector of 3 elements
        h11 = isi_channel(nu);
        h12 = isi_channel(nu);
%         h21 = isi_channel(nu);
%         h22 = isi_channel(nu);

        %this function creates a lower toeplitz square matrix of dimension (n+nu)
        %whose first column is the nu+1 impulse responses of h11
%         H11 = lower_toeplitz(n, nu, h11);
%         H12 = lower_toeplitz(n, nu, h12);
        H11 = circulant_matrix(n, 0, h11); % put nu=0 here as H11 is n by n for SC-FDE
        H12 = circulant_matrix(n, 0, h12); % put nu=0 here as H11 is n by n for SC-FDE
%         H21 = lower_toeplitz(n, nu, h21);
%         H22 = lower_toeplitz(n, nu, h22);
        
        
        lemda11 = Q * H11 * Q';
        lemda12 = Q * H12 * Q';
        
%         lemda11 = conj(Q) * H11 * Q;
%         lemda12 = conj(Q) * H12 * Q;
        
        lemda = [lemda11 lemda12 ; lemda12' -(lemda11')];


        %Noise generation
        noise1 = sqrt(1/snr) .* (randn(n, 1) + j * randn(n, 1)) ./ sqrt(2);
        noise2 = sqrt(1/snr) .* (randn(n, 1) + j * randn(n, 1)) ./ sqrt(2);
        
%         noise1 = 10 * exp(-10) * ones(n,1);
%         noise2 = 10 * exp(-10) * ones(n,1);

        y1 = H11 * C11 + H12 * C21 + noise1;
%         y2 = - H11 * (conj(C21)) + H12 * (conj(C11)) + noise2;
        y2 = - H11' * C21 + H12' * C11 + conj(noise2);

        Y1 = Q * y1;
        Y2 = Q * y2;
                
%         Y = [Y1 ; conj(Y2)];
        Y = [Y1 ; Y2];
        
        Y_tilda = lemda' * Y;
%         check1 = Y_tilda - lemda' * lemda * blkdiag(Q,Q) * [C11 ; C21]; %fine
        
        lem = lemda' * lemda;
        lemda_0 = lem(1:n,1:n);
%         check2 = Y_tilda(1:n) - lemda_0 * Q * C11; %fine
%         check3 = Y_tilda(n+1:2*n) - lemda_0 * Q * C21; % fine
% %         Z_zfe = [inv(lemda_0) * Y_tilda(1:n) ; inv(lemda_0) * Y_tilda(n+1:2*n)];
        Z_zfe = inv(lemda) * Y;
        
        z_zfe = [Q' * Z_zfe(1:n) ; Q' * Z_zfe(n+1:2*n)];
        
        Signal_at_SC_FDE = [lemda_0 * Q * C11 ; lemda_0 * Q * C21];
        Noise_at_SC_FDE = lemda' * [Q * noise1 ; conj(Q * noise2)];
        SNR_at_SC_FDE = (Signal_at_SC_FDE' * Signal_at_SC_FDE) / (Noise_at_SC_FDE' * Noise_at_SC_FDE);
        
        lemda_tilda = inv(lemda' * lemda + (1/SNR_at_SC_FDE) * eye(2*n));
        
        Z_mmse_le = lemda_tilda * Y_tilda;
        
        z_mmse_le = [Q' * Z_mmse_le(1:n) ; Q' * Z_mmse_le(n+1:2*n)];
        

        
        % ZF-DFE before separating the 2 sub-blocks gives the same
        % performance as after separating the 2 sub-blocks
        
        [P_zf_dfe,R_zf_dfe] = qr(lemda * blkdiag(Q,Q));
        
        y_match = (lemda * blkdiag(Q,Q))' * Y;
        
        Noise_at_mmse_dfe = (lemda * blkdiag(Q,Q))' * [Q * noise1 ; conj(Q * noise2)];
        Signal_at_mmse_dfe = y_match - Noise_at_mmse_dfe;
        SNR_at_mmse_dfe = (Signal_at_mmse_dfe' * Signal_at_mmse_dfe) / (Noise_at_mmse_dfe' * Noise_at_mmse_dfe);
        
        LEMDA_mmse = [lemda * blkdiag(Q,Q) ; sqrt(1/SNR_at_mmse_dfe) * eye(2*n)];
        [P_mmse_dfe, R_mmse_dfe] = qr(LEMDA_mmse);
        R_mmse_dfe = R_mmse_dfe(1:2*n,1:2*n);
        
        y_zf_dfe = P_zf_dfe' * Y;
        y_mmse_dfe = inv(R_mmse_dfe') * y_match;
        
        z_zf_dfe(2*n) = y_zf_dfe(2*n)/R_zf_dfe(2*n,2*n);
        z_zf_dfe_bits = [sign(real(z_zf_dfe(2*n)))+1 sign(imag(z_zf_dfe(2*n)))+1]./2;
        z_zf_dfe(2*n) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(2 * n);
        
        z_mmse_dfe(2*n) = y_mmse_dfe(2*n)/R_mmse_dfe(2*n,2*n);
        z_mmse_dfe_bits = [sign(real(z_mmse_dfe(2*n)))+1 sign(imag(z_mmse_dfe(2*n)))+1]./2;
        z_mmse_dfe(2*n) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(2 * n);
        
        for kk = 2*n-1:-1:1
            isi1_zf_dfe = 0;
            isi1_mmse_dfe = 0;

            for ll = kk + 1:2*n
                isi1_zf_dfe = isi1_zf_dfe + R_zf_dfe(kk,ll) * z_zf_dfe(ll);
                isi1_mmse_dfe = isi1_mmse_dfe + R_mmse_dfe(kk,ll) * z_mmse_dfe(ll);
            end

            z_zf_dfe(kk) = (y_zf_dfe(kk) - isi1_zf_dfe)/R_zf_dfe(kk,kk);
            z_zf_dfe_bits = [sign(real(z_zf_dfe(kk)))+1 sign(imag(z_zf_dfe(kk)))+1]./2;
            z_zf_dfe(kk) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(2 * n);
            
            z_mmse_dfe(kk) = (y_mmse_dfe(kk) - isi1_mmse_dfe)/R_mmse_dfe(kk,kk);
            z_mmse_dfe_bits = [sign(real(z_mmse_dfe(kk)))+1 sign(imag(z_mmse_dfe(kk)))+1]./2;
            z_mmse_dfe(kk) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(2 * n);
        end

        % %         % QR factorization
% %         
% %         [P,R] = qr(lemda_0 * Q);
% %         
% %         z1_z = P' * Y_tilda(1:n);
% %         z2_z = P' * Y_tilda(n+1:2*n);
% %         
% %         %MMSE-DFE
% %         
% %         LEMDA_mmse = [lemda_0 * Q ; sqrt(1/SNR_at_SC_FDE) * eye(n)];
% %         [P_mmse_dfe, R_mmse_dfe] = qr(LEMDA_mmse);
% % %         check4 = Y_tilda(1:n) - P_mmse_dfe(1:n,1:n) * R(1:n,:) * C11; %fine
% % %         check5 = Y_tilda(n+1:2*n) - P_mmse_dfe(1:n,1:n) * R(1:n,:) * C21; %fine
% %         
% %         R_mmse_dfe = R_mmse_dfe(1:n,:);
% %         
% %         % Now P_mmse_dfe is 2n by 2n and R_mmse_dfe is 2n by n.
% %         % P_mmse_dfe' * P_mmse_dfe = identity but
% %         % P_mmse_dfe(1:n,1:n)' * P_mmse_dfe(1:n,1:n) != identity
% %         z1_m = inv(P_mmse_dfe(1:n,1:n)) * Y_tilda(1:n);
% %         z2_m = inv(P_mmse_dfe(1:n,1:n)) * Y_tilda(n+1:2*n);

%         % Z_zf_dfe = D^(-1/2) * Z;
% % %         z1_zf_dfe(n) = z1_z(n)/R(n,n);
% % %         z2_zf_dfe(n) = z2_z(n)/R(n,n);
% % %         z1_mmse_dfe(n) = z1_m(n)/R_mmse_dfe(n,n);
% % %         z2_mmse_dfe(n) = z2_m(n)/R_mmse_dfe(n,n);
% % %         
% % %         z1_zf_dfe_bits = [sign(real(z1_zf_dfe(n)))+1 sign(imag(z1_zf_dfe(n)))+1]./2;
% % %         z2_zf_dfe_bits = [sign(real(z2_zf_dfe(n)))+1 sign(imag(z2_zf_dfe(n)))+1]./2;
% % %         z1_mmse_dfe_bits = [sign(real(z1_mmse_dfe(n)))+1 sign(imag(z1_mmse_dfe(n)))+1]./2;
% % %         z2_mmse_dfe_bits = [sign(real(z2_mmse_dfe(n)))+1 sign(imag(z2_mmse_dfe(n)))+1]./2;
% % % %         
% % % %         % All the symbols are divided by sqrt(2*n) to make them same as the
% % % %         % transmitted symbols C11 and C21
% % %         z1_zf_dfe(n) = qpsk_from_bits(z1_zf_dfe_bits)/sqrt(2 * n);
% % %         z2_zf_dfe(n) = qpsk_from_bits(z2_zf_dfe_bits)/sqrt(2 * n);
% % %         z1_mmse_dfe(n) = qpsk_from_bits(z1_mmse_dfe_bits)/sqrt(2 * n);
% % %         z2_mmse_dfe(n) = qpsk_from_bits(z2_mmse_dfe_bits)/sqrt(2 * n);
% % % % 
% % %         for kk = n-1:-1:1
% % %             isi1_zf_dfe = 0;
% % %             isi2_zf_dfe = 0;
% % %             isi1_mmse_dfe = 0;
% % %             isi2_mmse_dfe = 0;
% % % 
% % % %             memory = min(n - kk, nu);
% % % 
% % %             for ll = kk + 1:n
% % %                 isi1_zf_dfe = isi1_zf_dfe + R(kk,ll) * z1_zf_dfe(ll);
% % %                 isi2_zf_dfe = isi2_zf_dfe + R(kk,ll) * z2_zf_dfe(ll);
% % %                 isi1_mmse_dfe = isi1_mmse_dfe + R_mmse_dfe(kk,ll) * z1_mmse_dfe(ll);
% % %                 isi2_mmse_dfe = isi2_mmse_dfe + R_mmse_dfe(kk,ll) * z2_mmse_dfe(ll);
% % % % keyboard
% % %             end
% % % 
% % %             z1_zf_dfe(kk) = (z1_z(kk) - isi1_zf_dfe)/R(kk,kk);
% % %             z2_zf_dfe(kk) = (z2_z(kk) - isi2_zf_dfe)/R(kk,kk);
% % %             z1_mmse_dfe(kk) = (z1_m(kk) - isi1_mmse_dfe)/R_mmse_dfe(kk,kk);
% % %             z2_mmse_dfe(kk) = (z2_m(kk) - isi2_mmse_dfe)/R_mmse_dfe(kk,kk);
% % %               
% % %             z1_zf_dfe_bits = [sign(real(z1_zf_dfe(kk)))+1 sign(imag(z1_zf_dfe(kk)))+1]./2;
% % %             z2_zf_dfe_bits = [sign(real(z2_zf_dfe(kk)))+1 sign(imag(z2_zf_dfe(kk)))+1]./2;
% % %             z1_mmse_dfe_bits = [sign(real(z1_mmse_dfe(kk)))+1 sign(imag(z1_mmse_dfe(kk)))+1]./2;
% % %             z2_mmse_dfe_bits = [sign(real(z2_mmse_dfe(kk)))+1 sign(imag(z2_mmse_dfe(kk)))+1]./2;
% % %             
% % %             z1_zf_dfe(kk) = qpsk_from_bits(z1_zf_dfe_bits)/sqrt(2 * n);
% % %             z2_zf_dfe(kk) = qpsk_from_bits(z2_zf_dfe_bits)/sqrt(2 * n);
% % %             z1_mmse_dfe(kk) = qpsk_from_bits(z1_mmse_dfe_bits)/sqrt(2 * n);
% % %             z2_mmse_dfe(kk) = qpsk_from_bits(z2_mmse_dfe_bits)/sqrt(2 * n);
% % %         end
%         

%         
%         %MMSE-OSIC
        Signal_at_mmse_osic = Y - Noise_at_SC_FDE;
        SNR_at_mmse_osic = (Signal_at_mmse_osic' * Signal_at_mmse_osic) / (Noise_at_SC_FDE' * Noise_at_SC_FDE);
        z_mmse_osic = osic_mmse(lemda * blkdiag(Q,Q), Y, SNR_at_mmse_osic, 2*n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555

        %Decoding
  
        rec_bits_zfe = qpsk_decoded_bits(z_zfe); %([Z1_zfe ; Z2_zfe]);
        rec_bits_mmse_le = qpsk_decoded_bits(z_mmse_le);
        rec_bits_zf_dfe = qpsk_decoded_bits(z_zf_dfe.');
        rec_bits_mmse_dfe = qpsk_decoded_bits(z_mmse_dfe.');
        
        rec_bits_mmse_osic = qpsk_decoded_bits(z_mmse_osic);
        
        error_zfe = (rec_bits_zfe ~= trans_bits);
        error_mmse_le = (rec_bits_mmse_le ~= trans_bits);
        error_zf_dfe = (rec_bits_zf_dfe ~= trans_bits);
        error_mmse_dfe = (rec_bits_mmse_dfe ~= trans_bits);
        error_mmse_osic = (rec_bits_mmse_osic ~= trans_bits);

        total_error_zfe = total_error_zfe + sum(error_zfe);
        total_error_mmse_le = total_error_mmse_le + sum(error_mmse_le);
        total_error_zf_dfe = total_error_zf_dfe + sum(error_zf_dfe);
        total_error_mmse_dfe = total_error_mmse_dfe + sum(error_mmse_dfe);
        total_error_mmse_osic = total_error_mmse_osic + sum(error_mmse_osic);
        
        if mod(l,20000)==0
%             k
            l
        end

    end

    BER_zfe(k) = total_error_zfe / (2 * n * bits_per_symbol * run);
    BER_mmse_le(k) = total_error_mmse_le /(2 * n * bits_per_symbol * run);
    BER_zf_dfe(k) = total_error_zf_dfe /(2 * n * bits_per_symbol * run);
    BER_mmse_dfe(k) = total_error_mmse_dfe /(2 * n * bits_per_symbol * run);
    BER_mmse_osic(k) = total_error_mmse_osic /(2 * n * bits_per_symbol * run);
end



% semilogy(SNR,BER_zfe,'b*--',SNR,BER_mmse_le,'ms--',SNR,BER_zf_dfe,'ro:',SNR, BER_mmse_dfe,'y^--', SNR, BER_mmse_osic, 'gv-.')
semilogy(SNR,BER_zfe,'b*:',SNR,BER_mmse_le,'ms--', SNR,BER_zf_dfe,'ro:', SNR, BER_mmse_dfe,'gv:',SNR, BER_mmse_osic,'kx:')
legend('ZFE', 'MMSE-LE', 'ZF-DFE','MMSE-DFE','MMSE-OSIC');%, 'northeast');
grid;
% legend('ZFE', 'MMSE-LE', 'ZF-DFE', 'MMSE-DFE');%, 'northeast');
xlabel('SNR in dB');
ylabel('Bit Error Rate');
title('Performance comparison of different equalization techniques N = 4, nu = 1');
% title('Performance comparison of different equalization techniques for STBC in ISI multiple access channels');
% axis([5,25,10^(-6),10^0]);

