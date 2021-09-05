close all; clear all; clc;

% input bit energy
Ex = 1;

%input SNR in dB
SNR = [5:5:25];
snr_linear = 10.^(SNR/10);

num_of_user = 3;
num_of_Rx = 3;

% sub-block length is n i.e. one sub-block contains n symbols (symbols can
% be real as PAM or complex as QPSK
n = 4;
bits_per_symbol = 2;

%channel memory nu so channel will have nu+1 taps
nu = 1;

%number of transmit antenna
Mt = 2;

%changes the state of the random generator
randn('state', sum(100*clock));
rand('state', sum(100*clock));

%run = number of transmission
run = 100000;

%this matrix is called zero-stuffing matrix with dimension (n)*n.
Izs = [eye(n); zeros(nu,n)];

% The n by n DFT matrix
Q = DFT_matrix(n);

for k = 1:length(snr_linear)
    snr = snr_linear(k)
    
    total_error_zfe = zeros(1,num_of_user);
    total_error_mmse_le = zeros(1,num_of_user);
    total_error_zf_dfe = zeros(1,num_of_user);
    total_error_mmse_dfe = zeros(1,num_of_user);
    total_error_mmse_osic = zeros(1,num_of_user);
    
    total_error_zf_dfe_u = zeros(1,num_of_user);
    total_error_mmse_dfe_u = zeros(1,num_of_user);
    total_error_mmse_osic_u2 = zeros(1,num_of_user);
    
%     total_error_zf_dfe_u1 = 0;
%     total_error_zf_dfe_u2 = 0;
%     total_error_zf_dfe_u3 = 0;
    
      
%     total_error_mmse_osic = 0;

%     H1 = [];
%     H2 = [];
%     H3 = [];
    S = [];
    trans_bits = [];
    r = [];

    for l = 1:run
        H = [];
        noise = [];
        
        for u = 1:num_of_user
    
            % For 2 sub-blocks there are 2n symbols
            trans_bits{u} = randint(1,bits_per_symbol * Mt * n);
            sym = qpsk_from_bits(trans_bits{u});
            symbols = sym/sqrt(Mt * n);
            C11 = symbols(1:n).';
            C21 = symbols(n+1:2*n).';
            s{u} = [C11 ; C21];
            S{u} = [Q * C11 ; Q * C21];
            
            H_temp = [];
            
            for ii = 1:num_of_Rx
                for jj = 1:Mt
                    h = isi_channel(nu);
                    Hc{jj} = circulant_matrix(n,0,h);
                end
                H_temp = [H_temp ; Hc{1} Hc{2} ; -Hc{2}' Hc{1}'];
                
                %Noise generation. 1 and 2 for Rx 1 and 2 and 3 for Rx 2
                noise1 = sqrt(1/snr) .* (randn(n, 1) + j * randn(n, 1)) ./ sqrt(2);
                noise2 = sqrt(1/snr) .* (randn(n, 1) + j * randn(n, 1)) ./ sqrt(2);
            end
            
            H = [H H_temp];
            
            noise = [noise ; noise1 ; conj(noise2)];
            
% % % %             %ISI CHANNEL    
% % % % 
% % % %             % hkl specifies the channel from Tx antenna l to Rx antenna k
% % % %             % all are IID gaussian 3 tap channel so a vector of 3 elements
% % % %             h11 = isi_channel(nu);
% % % %             h12 = isi_channel(nu);
% % % %             h21 = isi_channel(nu);
% % % %             h22 = isi_channel(nu);
% % % %             h31 = isi_channel(nu);
% % % %             h32 = isi_channel(nu);
% % % %             
% % % %             % n dimentional square circulant matrix
% % % %             H11c = circulant_matrix(n, 0, h11);
% % % %             H12c = circulant_matrix(n, 0, h12);
% % % %             H21c = circulant_matrix(n, 0, h21);
% % % %             H22c = circulant_matrix(n, 0, h22);
% % % %             H31c = circulant_matrix(n, 0, h31);
% % % %             H32c = circulant_matrix(n, 0, h32);
% % % %             
% % % %             % H4(u) = overall channel matrix for user u at Rx antenna 4
% % % %             % assuming each user has 2 Tx antenna. 2*(n) dimensional
% % % %             % square matrix
% % % %             H1{u} = [H11c H12c ; H12c' -H11c'];
% % % %             H2{u} = [H21c H22c ; H22c' -H21c'];
% % % %             H3{u} = [H31c H32c ; H32c' -H31c'];
% % % %             
% % % % %             lemda_11{u} = Q * H11c * Q';
% % % % %             lemda_12{u} = Q * H12c * Q';
% % % % %             lemda_21{u} = Q * H21c * Q';
% % % % %             lemda_22{u} = Q * H22c * Q';
% % % % %             lemda_31{u} = Q * H31c * Q';
% % % % %             lemda_32{u} = Q * H32c * Q';
% % % % %             
% % % % %             lemda_111 = [lemda_11{u} lemda_12{u} ; lemda_12{u}' -lemda_11{u}'];
% % % % %             lemda_211 = [lemda_21{u} lemda_22{u} ; lemda_22{u}' -lemda_21{u}'];
% % % % %             lemda_311 = [lemda_31{u} lemda_32{u} ; lemda_32{u}' -lemda_31{u}'];
% % % % %             
% % % % %             lemda{u} = [lemda_111 ; lemda_211 ; lemda_311];
        end
            
%         lemda5 = [lemda{1} lemda{2} lemda{3}];
        
% % % %         %Noise generation. 1 and 2 for Rx 1 and 2 and 3 for Rx 2
% % % %         noise1 = sqrt(1/snr) .* (randn(n, 1) + j * randn(n, 1)) ./ sqrt(2);
% % % %         noise2 = sqrt(1/snr) .* (randn(n, 1) + j * randn(n, 1)) ./ sqrt(2);
% % % %         noise3 = sqrt(1/snr) .* (randn(n, 1) + j * randn(n, 1)) ./ sqrt(2);
% % % %         noise4 = sqrt(1/snr) .* (randn(n, 1) + j * randn(n, 1)) ./ sqrt(2);
% % % %         noise5 = sqrt(1/snr) .* (randn(n, 1) + j * randn(n, 1)) ./ sqrt(2);
% % % %         noise6 = sqrt(1/snr) .* (randn(n, 1) + j * randn(n, 1)) ./ sqrt(2);
        
% %         noise = [noise1 ; conj(noise2) ; noise3 ; conj(noise4) ; noise5 ; conj(noise6)];
% %         NOISE = blkdiag(Q,Q,Q,Q,Q,Q) * noise;
        NOISE = kron(eye(Mt * num_of_Rx),Q) * noise; % as num_of_Rx = num_of_user
        
        % Overall channnel matrix
% %         H = [H1{1} H1{2} H1{3} ; H2{1} H2{2} H2{3} ; H3{1} H3{2} H3{3}];
% %         LEMDA = blkdiag(Q,Q,Q,Q,Q,Q) * H * blkdiag(Q',Q',Q',Q',Q',Q');
        LEMDA = kron(eye(Mt * num_of_Rx),Q) * H * kron(eye(Mt * num_of_Rx),Q');

        s_tilda = [s{1} ; s{2} ; s{3}];
        
        % Total received signal at the receiver
        y_total = H * s_tilda + noise;
%         Y_total = blkdiag(Q,Q,Q,Q,Q,Q) * y_total;
        Y_total = kron(eye(Mt * num_of_Rx),Q) * y_total;
        Signal_at_mmse_osic = Y_total - NOISE;
        SNR_at_mmse_osic = (Signal_at_mmse_osic' * Signal_at_mmse_osic) / (NOISE' * NOISE);

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ZF-DFE before decoupling the users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        for u = 1:num_of_user
            
            % if u = 1 no shift is reqd, after that we need a shift of 2n
            if u == 1
                LEMDA_u = LEMDA;
            else
                LEMDA_u = circshift(LEMDA_u,[0, -2*n]);
            end
            
            
%         [P_zf_dfe_u,R_zf_dfe_u] = qr(LEMDA * blkdiag(Q,Q,Q,Q,Q,Q));
            [P_zf_dfe_u,R_zf_dfe_u] = qr(LEMDA_u * kron(eye(Mt * num_of_Rx),Q));
        
            % MMSE-DFE
            y_mmse_match = (LEMDA_u * kron(eye(Mt * num_of_Rx),Q))' * Y_total;
            
            Noise_at_mmse_dfe = (LEMDA_u * kron(eye(Mt * num_of_Rx),Q))' * NOISE;
            Signal_at_mmse_dfe = y_mmse_match - Noise_at_mmse_dfe;
            SNR_at_mmse_dfe = (Signal_at_mmse_dfe' * Signal_at_mmse_dfe) / (Noise_at_mmse_dfe' * Noise_at_mmse_dfe);
           
            [P_mmse_dfe_u, R_mmse_dfe_u] = qr([LEMDA_u * kron(eye(Mt * num_of_Rx),Q); sqrt(1/SNR_at_mmse_dfe) * eye(2*num_of_user*n)]);
            R_mmse_dfe_u = R_mmse_dfe_u(1:2*num_of_user*n,:);

            y_zf_dfe_u = P_zf_dfe_u' * Y_total;
            y_mmse_dfe_u = inv(R_mmse_dfe_u') * y_mmse_match;
        
            z_zf_dfe_u(2*num_of_user*n) = y_zf_dfe_u(2*num_of_user*n)/R_zf_dfe_u(2*num_of_user*n,2*num_of_user*n);
            z_zf_dfe_bits = [sign(real(z_zf_dfe_u(2*num_of_user*n)))+1 sign(imag(z_zf_dfe_u(2*num_of_user*n)))+1]./2;
            z_zf_dfe_u(2*num_of_user*n) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(Mt * n);
        
            z_mmse_dfe_u(2*num_of_user*n) = y_mmse_dfe_u(2*num_of_user*n)/R_mmse_dfe_u(2*num_of_user*n,2*num_of_user*n);
            z_mmse_dfe_bits = [sign(real(z_mmse_dfe_u(2*num_of_user*n)))+1 sign(imag(z_mmse_dfe_u(2*num_of_user*n)))+1]./2;
            z_mmse_dfe_u(2*num_of_user*n) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(Mt * n);
        
            for kk = 2*num_of_user*n-1:-1:1
                isi1_zf_dfe_u = 0;
                isi1_mmse_dfe_u = 0;
%               isi2_zf_dfe_u = 0;
%               isi1_mmse_dfe = 0;
%               isi2_mmse_dfe = 0;

                for ll = kk + 1:2*num_of_user*n
                    isi1_zf_dfe_u = isi1_zf_dfe_u + R_zf_dfe_u(kk,ll) * z_zf_dfe_u(ll);
%                 isi2_zf_dfe_u = isi2_zf_dfe_u + R2(kk,ll) * z_zf_dfe_u_2(ll);
                    isi1_mmse_dfe_u = isi1_mmse_dfe_u + R_mmse_dfe_u(kk,ll) * z_mmse_dfe_u(ll);
%                 isi1_mmse_dfe = isi1_mmse_dfe + R_mmse_dfe(kk,ll) * z_mmse_dfe(ll);
%                 isi2_mmse_dfe = isi2_mmse_dfe + R_mmse_dfe_2(kk,ll) * z_mmse_dfe_2(ll);
                end

                z_zf_dfe_u(kk) = (y_zf_dfe_u(kk) - isi1_zf_dfe_u)/R_zf_dfe_u(kk,kk);
                z_zf_dfe_bits = [sign(real(z_zf_dfe_u(kk)))+1 sign(imag(z_zf_dfe_u(kk)))+1]./2;
                z_zf_dfe_u(kk) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(Mt * n);

%               z_zf_dfe_u_2(kk) = (y_zf_dfe_u_2(kk) - isi2_zf_dfe_u)/R2(kk,kk);
%               z_zf_dfe_bits = [sign(real(z_zf_dfe_u_2(kk)))+1 sign(imag(z_zf_dfe_u_2(kk)))+1]./2;
%               z_zf_dfe_u_2(kk) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(2 * n);
%             
                z_mmse_dfe_u(kk) = (y_mmse_dfe_u(kk) - isi1_mmse_dfe_u)/R_mmse_dfe_u(kk,kk);
                z_mmse_dfe_bits = [sign(real(z_mmse_dfe_u(kk)))+1 sign(imag(z_mmse_dfe_u(kk)))+1]./2;
                z_mmse_dfe_u(kk) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(2 * n);
%             
%               z_mmse_dfe_2(kk) = (y_mmse_dfe_2(kk) - isi2_mmse_dfe)/R_mmse_dfe_2(kk,kk);
%               z_mmse_dfe_bits = [sign(real(z_mmse_dfe_2(kk)))+1 sign(imag(z_mmse_dfe_2(kk)))+1]./2;
%               z_mmse_dfe_2(kk) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(2 * n);
            end
            
            z_mmse_osic_u2 = osic_mmse(LEMDA_u * kron(eye(Mt * num_of_Rx),Q), Y_total,...
                             SNR_at_mmse_osic, Mt*n);
        
            rec_bits_zf_dfe_u = qpsk_decoded_bits(z_zf_dfe_u(1:2*n).');
            rec_bits_mmse_dfe_u = qpsk_decoded_bits(z_mmse_dfe_u(1:2*n).');
            rec_bits_mmse_osic_u2 = qpsk_decoded_bits(z_mmse_osic_u2(4*n+1:6*n));
            
            error_zf_dfe_u = (rec_bits_zf_dfe_u ~= trans_bits{u});
            error_mmse_dfe_u = (rec_bits_mmse_dfe_u ~= trans_bits{u});
            
            if u == 1
                error_mmse_osic_u2 = (rec_bits_mmse_osic_u2 ~= trans_bits{num_of_user});
                total_error_mmse_osic_u2(num_of_user) = total_error_mmse_osic_u2(num_of_user) + sum(error_mmse_osic_u2);
            else
                error_mmse_osic_u2 = (rec_bits_mmse_osic_u2 ~= trans_bits{u-1});
                total_error_mmse_osic_u2(u-1) = total_error_mmse_osic_u2(u-1) + sum(error_mmse_osic_u2);
            end
            
            total_error_zf_dfe_u(u) = total_error_zf_dfe_u(u) + sum(error_zf_dfe_u);
            total_error_mmse_dfe_u(u) = total_error_mmse_dfe_u(u) + sum(error_mmse_dfe_u);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of DFE before decoupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for uu = num_of_user:-1:2

            A = LEMDA(1:2*(uu-1)*(n),1:2*(uu-1)*(n));
            B = LEMDA(1:2*(uu-1)*(n),2*(uu-1)*(n)+1:2*uu*(n));
            C = LEMDA(2*(uu-1)*(n)+1:2*(uu-1)*(n)+2*(n),1:2*(uu-1)*(n));
            D = LEMDA(2*(uu-1)*(n)+1:2*(uu-1)*(n)+2*(n),2*(uu-1)*(n)+1:2*uu*(n));


            W_decorr = [eye(2*(uu-1)*(n)) -B * inv(D) ; -C * inv(A) eye(2*(n))];
%             temp = W_decorr * LEMDA;
            Y_total = W_decorr * Y_total;
            Y{uu} = Y_total(2*(uu-1)*(n)+1:2*uu*(n));
            Y_total = Y_total(1:2*(uu-1)*(n));

            N_decorr = W_decorr * NOISE;
            NOISE = N_decorr(1:2*(uu-1)*n);

            LEMDA = A - B * inv(D) * C;
            DELTA = D - C * inv(A) * B;

            % Matched filter output
            Y_match = DELTA' * Y{uu};

            Noise_at_SC_FDE = DELTA' * N_decorr(2*(uu-1)*(n)+1:2*uu*(n));
            Signal_at_SC_FDE = Y_match - Noise_at_SC_FDE;
            SNR_at_SC_FDE = (Signal_at_SC_FDE' * Signal_at_SC_FDE) / (Noise_at_SC_FDE' * Noise_at_SC_FDE);

            % MMSE-LE
            W_mmse_le = inv(DELTA' * DELTA + (1/SNR_at_SC_FDE) * eye(2*n));

            % Output of ZFE
            Z_zfe = inv(DELTA) * Y{uu};
            z_zfe = kron(eye(Mt),Q') * [Z_zfe(1:n) ; Z_zfe(n+1:2*n)];

            % Output of MMSE-LE
            Z_mmse_le = W_mmse_le * Y_match;
            z_mmse_le = kron(eye(Mt),Q') * [Z_mmse_le(1:n) ; Z_mmse_le(n+1:2*n)];
            
            %MMSE-OSIC
            z_mmse_osic = osic_mmse(DELTA * blkdiag(Q,Q), Y{uu}, SNR_at_SC_FDE, 2*n);

            % ZF-DFE
            [P_zf_dfe,R_zf_dfe] = qr(DELTA * blkdiag(Q,Q));

            y_zf_dfe = P_zf_dfe' * Y{uu};

            z_zf_dfe(2*n) = y_zf_dfe(2*n)/R_zf_dfe(2*n,2*n);
            z_zf_dfe_bits = [sign(real(z_zf_dfe(2*n)))+1 sign(imag(z_zf_dfe(2*n)))+1]./2;
            z_zf_dfe(2*n) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(Mt * n);

            % MMSE-DFE
            y_mmse_match = (DELTA * blkdiag(Q,Q))' * Y{uu};

            Noise_at_mmse_dfe = (DELTA * blkdiag(Q,Q))' * N_decorr(2*(uu-1)*(n)+1:2*uu*(n));
            Signal_at_mmse_dfe = y_mmse_match - Noise_at_mmse_dfe;
            SNR_at_mmse_dfe = (Signal_at_mmse_dfe' * Signal_at_mmse_dfe) / (Noise_at_mmse_dfe' * Noise_at_mmse_dfe);

            LEMDA_mmse = [DELTA * blkdiag(Q,Q) ; sqrt(1/SNR_at_mmse_dfe) * eye(2*n)];
            [P_mmse_dfe, R_mmse_dfe] = qr(LEMDA_mmse);
            R_mmse_dfe = R_mmse_dfe(1:2*n,1:2*n);

            y_mmse_dfe = inv(R_mmse_dfe') * y_mmse_match;

            z_mmse_dfe(2*n) = y_mmse_dfe(2*n)/R_mmse_dfe(2*n,2*n);
            z_mmse_dfe_bits = [sign(real(z_mmse_dfe(2*n)))+1 sign(imag(z_mmse_dfe(2*n)))+1]./2;
            z_mmse_dfe(2*n) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(Mt * n);

            for kk = 2*n-1:-1:1
            isi1_zf_dfe = 0;
%             isi2_zf_dfe = 0;
            isi1_mmse_dfe = 0;
%             isi2_mmse_dfe = 0;

            for ll = kk + 1:2*n
                isi1_zf_dfe = isi1_zf_dfe + R_zf_dfe(kk,ll) * z_zf_dfe(ll);
%                 isi2_zf_dfe = isi2_zf_dfe + R2(kk,ll) * z_zf_dfe_2(ll);

                isi1_mmse_dfe = isi1_mmse_dfe + R_mmse_dfe(kk,ll) * z_mmse_dfe(ll);
%                 isi2_mmse_dfe = isi2_mmse_dfe + R_mmse_dfe_2(kk,ll) * z_mmse_dfe_2(ll);
            end

            z_zf_dfe(kk) = (y_zf_dfe(kk) - isi1_zf_dfe)/R_zf_dfe(kk,kk);
            z_zf_dfe_bits = [sign(real(z_zf_dfe(kk)))+1 sign(imag(z_zf_dfe(kk)))+1]./2;
            z_zf_dfe(kk) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(Mt * n);

            z_mmse_dfe(kk) = (y_mmse_dfe(kk) - isi1_mmse_dfe)/R_mmse_dfe(kk,kk);
            z_mmse_dfe_bits = [sign(real(z_mmse_dfe(kk)))+1 sign(imag(z_mmse_dfe(kk)))+1]./2;
            z_mmse_dfe(kk) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(Mt * n);

%             z_zf_dfe_2(kk) = (y_zf_dfe_2(kk) - isi2_zf_dfe)/R2(kk,kk);
%             z_zf_dfe_bits = [sign(real(z_zf_dfe_2(kk)))+1 sign(imag(z_zf_dfe_2(kk)))+1]./2;
%             z_zf_dfe_2(kk) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(2 * n);
%             
%             z_mmse_dfe(kk) = (y_mmse_dfe(kk) - isi1_mmse_dfe)/R_mmse_dfe(kk,kk);
%             z_mmse_dfe_bits = [sign(real(z_mmse_dfe(kk)))+1 sign(imag(z_mmse_dfe(kk)))+1]./2;
%             z_mmse_dfe(kk) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(2 * n);
%             
%             z_mmse_dfe_2(kk) = (y_mmse_dfe_2(kk) - isi2_mmse_dfe)/R_mmse_dfe_2(kk,kk);
%             z_mmse_dfe_bits = [sign(real(z_mmse_dfe_2(kk)))+1 sign(imag(z_mmse_dfe_2(kk)))+1]./2;
%             z_mmse_dfe_2(kk) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(2 * n);
            end

            % Decoding            
            rec_bits_zfe = qpsk_decoded_bits(z_zfe);
            rec_bits_mmse_le = qpsk_decoded_bits(z_mmse_le);
            rec_bits_zf_dfe = qpsk_decoded_bits(z_zf_dfe.');
            rec_bits_mmse_dfe = qpsk_decoded_bits(z_mmse_dfe.');
            rec_bits_mmse_osic = qpsk_decoded_bits(z_mmse_osic);

            error_zfe = (rec_bits_zfe ~= trans_bits{uu});
            error_mmse_le = (rec_bits_mmse_le ~= trans_bits{uu});
            error_zf_dfe = (rec_bits_zf_dfe ~= trans_bits{uu});
            error_mmse_dfe = (rec_bits_mmse_dfe ~= trans_bits{uu});
            error_mmse_osic = (rec_bits_mmse_osic ~= trans_bits{uu});

            total_error_zfe(uu) = total_error_zfe(uu) + sum(error_zfe);
            total_error_mmse_le(uu) = total_error_mmse_le(uu) + sum(error_mmse_le);
            total_error_zf_dfe(uu) = total_error_zf_dfe(uu) + sum(error_zf_dfe);
            total_error_mmse_dfe(uu) = total_error_mmse_dfe(uu) + sum(error_mmse_dfe);
            total_error_mmse_osic(uu) = total_error_mmse_osic(uu) + sum(error_mmse_osic);
            
        end
                
        
%%%%%%%%%%%%%%% Now process the first user

        Y_match = LEMDA' * Y_total;
        
        Noise_at_SC_FDE = LEMDA' * NOISE;
        Signal_at_SC_FDE = Y_match - Noise_at_SC_FDE;
        SNR_at_SC_FDE = (Signal_at_SC_FDE' * Signal_at_SC_FDE) / (Noise_at_SC_FDE' * Noise_at_SC_FDE);

        W_mmse_le = inv(LEMDA' * LEMDA + (1/SNR_at_SC_FDE) * eye(2*n));
        
        Z_zfe = inv(LEMDA) * Y_total;
        z_zfe = [Q' * Z_zfe(1:n) ; Q' * Z_zfe(n+1:2*n)];
        
        Z_mmse_le = W_mmse_le * Y_match;
        z_mmse_le = [Q' * Z_mmse_le(1:n) ; Q' * Z_mmse_le(n+1:2*n)];
        
        z_mmse_osic = osic_mmse(LEMDA * blkdiag(Q,Q), Y_total, SNR_at_SC_FDE, 2*n);
        
        % ZF-DFE
        [P_zf_dfe,R_zf_dfe] = qr(LEMDA * blkdiag(Q,Q));
        
        % MMSE-DFE
        y_mmse_match = (LEMDA * blkdiag(Q,Q))' * Y_total;
            
        Noise_at_mmse_dfe = (LEMDA * blkdiag(Q,Q))' * NOISE;
        Signal_at_mmse_dfe = y_mmse_match - Noise_at_mmse_dfe;
        SNR_at_mmse_dfe = (Signal_at_mmse_dfe' * Signal_at_mmse_dfe) / (Noise_at_mmse_dfe' * Noise_at_mmse_dfe);
            
        LEMDA_mmse = [LEMDA * blkdiag(Q,Q) ; sqrt(1/SNR_at_mmse_dfe) * eye(2*n)];
        [P_mmse_dfe, R_mmse_dfe] = qr(LEMDA_mmse);
        R_mmse_dfe = R_mmse_dfe(1:2*n,1:2*n);
            
        y_zf_dfe = P_zf_dfe' * Y_total;
        y_mmse_dfe = inv(R_mmse_dfe') * y_mmse_match;

        z_zf_dfe(2*n) = y_zf_dfe(2*n)/R_zf_dfe(2*n,2*n);
        z_zf_dfe_bits = [sign(real(z_zf_dfe(2*n)))+1 sign(imag(z_zf_dfe(2*n)))+1]./2;
        z_zf_dfe(2*n) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(Mt * n);
        
        z_mmse_dfe(2*n) = y_mmse_dfe(2*n)/R_mmse_dfe(2*n,2*n);
        z_mmse_dfe_bits = [sign(real(z_mmse_dfe(2*n)))+1 sign(imag(z_mmse_dfe(2*n)))+1]./2;
        z_mmse_dfe(2*n) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(Mt * n);
        
        for kk = 2*n-1:-1:1
            isi1_zf_dfe = 0;
%             isi2_zf_dfe = 0;
            isi1_mmse_dfe = 0;
%             isi2_mmse_dfe = 0;

            for ll = kk + 1:2*n
                isi1_zf_dfe = isi1_zf_dfe + R_zf_dfe(kk,ll) * z_zf_dfe(ll);
%                 isi2_zf_dfe = isi2_zf_dfe + R2(kk,ll) * z_zf_dfe_2(ll);
                
                isi1_mmse_dfe = isi1_mmse_dfe + R_mmse_dfe(kk,ll) * z_mmse_dfe(ll);
%                 isi2_mmse_dfe = isi2_mmse_dfe + R_mmse_dfe_2(kk,ll) * z_mmse_dfe_2(ll);
            end

            z_zf_dfe(kk) = (y_zf_dfe(kk) - isi1_zf_dfe)/R_zf_dfe(kk,kk);
            z_zf_dfe_bits = [sign(real(z_zf_dfe(kk)))+1 sign(imag(z_zf_dfe(kk)))+1]./2;
            z_zf_dfe(kk) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(Mt * n);
            
            z_mmse_dfe(kk) = (y_mmse_dfe(kk) - isi1_mmse_dfe)/R_mmse_dfe(kk,kk);
            z_mmse_dfe_bits = [sign(real(z_mmse_dfe(kk)))+1 sign(imag(z_mmse_dfe(kk)))+1]./2;
            z_mmse_dfe(kk) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(Mt * n);
            
%             z_zf_dfe_2(kk) = (y_zf_dfe_2(kk) - isi2_zf_dfe)/R2(kk,kk);
%             z_zf_dfe_bits = [sign(real(z_zf_dfe_2(kk)))+1 sign(imag(z_zf_dfe_2(kk)))+1]./2;
%             z_zf_dfe_2(kk) = qpsk_from_bits(z_zf_dfe_bits)/sqrt(2 * n);
%             
%             z_mmse_dfe(kk) = (y_mmse_dfe(kk) - isi1_mmse_dfe)/R_mmse_dfe(kk,kk);
%             z_mmse_dfe_bits = [sign(real(z_mmse_dfe(kk)))+1 sign(imag(z_mmse_dfe(kk)))+1]./2;
%             z_mmse_dfe(kk) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(2 * n);
%             
%             z_mmse_dfe_2(kk) = (y_mmse_dfe_2(kk) - isi2_mmse_dfe)/R_mmse_dfe_2(kk,kk);
%             z_mmse_dfe_bits = [sign(real(z_mmse_dfe_2(kk)))+1 sign(imag(z_mmse_dfe_2(kk)))+1]./2;
%             z_mmse_dfe_2(kk) = qpsk_from_bits(z_mmse_dfe_bits)/sqrt(2 * n);
            end


        % Decoding            
        rec_bits_zfe = qpsk_decoded_bits(z_zfe);
        rec_bits_mmse_le = qpsk_decoded_bits(z_mmse_le);
        rec_bits_zf_dfe = qpsk_decoded_bits(z_zf_dfe.');
        rec_bits_mmse_dfe = qpsk_decoded_bits(z_mmse_dfe.');
        rec_bits_mmse_osic = qpsk_decoded_bits(z_mmse_osic);

        error_zfe = (rec_bits_zfe ~= trans_bits{1});
        error_mmse_le = (rec_bits_mmse_le ~= trans_bits{1});
        error_zf_dfe = (rec_bits_zf_dfe ~= trans_bits{1});
        error_mmse_dfe = (rec_bits_mmse_dfe ~= trans_bits{1});
        error_mmse_osic = (rec_bits_mmse_osic ~= trans_bits{1});
        
% %         error_zf_dfe_u1 = (rec_bits_zf_dfe_u(1:2*bits_per_symbol*n) ~= trans_bits{1});
% %         error_zf_dfe_u2 = (rec_bits_zf_dfe_u(2*bits_per_symbol*n+1:4*bits_per_symbol*n) ~= trans_bits{2});
% %         error_zf_dfe_u3 = (rec_bits_zf_dfe_u(4*bits_per_symbol*n+1:6*bits_per_symbol*n) ~= trans_bits{3});
        
        total_error_zfe(uu-1) = total_error_zfe(uu-1) + sum(error_zfe);
        total_error_mmse_le(uu-1) = total_error_mmse_le(uu-1) + sum(error_mmse_le);
        total_error_zf_dfe(uu-1) = total_error_zf_dfe(uu-1) + sum(error_zf_dfe);
        total_error_mmse_dfe(uu-1) = total_error_mmse_dfe(uu-1) + sum(error_mmse_dfe);
        total_error_mmse_osic(uu-1) = total_error_mmse_osic(uu-1) + sum(error_mmse_osic);
        
%         total_error_zf_dfe_u1 = total_error_zf_dfe_u1 + sum(error_zf_dfe_u1);
%         total_error_zf_dfe_u2 = total_error_zf_dfe_u2 + sum(error_zf_dfe_u2);
%         total_error_zf_dfe_u3 = total_error_zf_dfe_u3 + sum(error_zf_dfe_u3);
        

        if mod(l,10000)==0
%             k
            l
        end
        
        
        
    end
    
    for uu = 1:num_of_user
        
        BER_zfe(uu,k) = total_error_zfe(uu) / (2 * n * bits_per_symbol * run);
        BER_mmse_le(uu,k) = total_error_mmse_le(uu) /(2 * n * bits_per_symbol * run);
        BER_zf_dfe(uu,k) = total_error_zf_dfe(uu) /(2 * n * bits_per_symbol * run);
        BER_mmse_dfe(uu,k) = total_error_mmse_dfe(uu) /(2 * n * bits_per_symbol * run);
        BER_mmse_osic(uu,k) = total_error_mmse_osic(uu) /(2 * n * bits_per_symbol * run);
        
        BER_zf_dfe_u(uu,k) = total_error_zf_dfe_u(uu) /(2 * n * bits_per_symbol * run);
        BER_mmse_dfe_u(uu,k) = total_error_mmse_dfe_u(uu) /(2 * n * bits_per_symbol * run);
        BER_mmse_osic_u2(uu,k) = total_error_mmse_osic_u2(uu) /(2 * n * bits_per_symbol * run);
        
    end
    
%     BER_zf_dfe_u1(k) = total_error_zf_dfe_u1 /(2 * n * bits_per_symbol * run);
        
%     BER_zfe_2(k) = total_error_zfe(2) / (2 * n * bits_per_symbol * run);
%     BER_mmse_le_2(k) = total_error_mmse_le(2) /(2 * n * bits_per_symbol * run);
%     BER_zf_dfe_2(k) = total_error_zf_dfe(2) /(2 * n * bits_per_symbol * run);
%     BER_mmse_dfe_2(k) = total_error_mmse_dfe(2) /(2 * n * bits_per_symbol * run);
    
%     BER_zf_dfe_u2(k) = total_error_zf_dfe_u2 /(2 * n * bits_per_symbol * run);
    
%     BER_zfe_3(k) = total_error_zfe(3) / (2 * n * bits_per_symbol * run);
%     BER_mmse_le_3(k) = total_error_mmse_le(3) /(2 * n * bits_per_symbol * run);
%     BER_zf_dfe_3(k) = total_error_zf_dfe(3) /(2 * n * bits_per_symbol * run);
%     BER_mmse_dfe_3(k) = total_error_mmse_dfe(3) /(2 * n * bits_per_symbol * run);
    
%     BER_zf_dfe_u3(k) = total_error_zf_dfe_u3 /(2 * n * bits_per_symbol * run);
    
%     BER_mmse_osic(k) = total_error_mmse_osic /(2 * n * bits_per_symbol * run);
end


% semilogy(SNR,BER_zfe_1,'b*--',SNR,BER_mmse_le_1,'ms--',SNR,BER_zfe_2,'ro:',SNR,BER_mmse_le_2,'gv:')
% semilogy(SNR,BER_zfe,'b*:',SNR,BER_mmse_le,'ms--',SNR, BER_mmse_osic, 'gv-.')
% semilogy(SNR,BER_zfe,'b*--',SNR,BER_mmse_le,'ms--',SNR,BER_zf_dfe,'ro:',SNR, BER_mmse_dfe,'y^--', SNR, BER_mmse_osic, 'gv-.')

for uu = 1:num_of_user
    
    figure(uu);
    set(gcf, 'Name','User'); % or set(figure(1), 'Name', 'User 1');
    % figure(aaa);
    % aaa = figure('Name', 'User1');
    % semilogy(SNR,BER_zfe_1,'b*--',SNR,BER_mmse_le_1,'ms--',SNR,BER_zf_dfe_1,'ro:',SNR, BER_mmse_dfe_1,'y^--')
    semilogy(SNR,BER_zfe(uu,:),'b*--',SNR,BER_mmse_le(uu,:),'ms--',SNR,BER_zf_dfe(uu,:),'ro:',...
             SNR, BER_mmse_dfe(uu,:),'g^--',SNR,BER_mmse_osic(uu,:),'k<:',SNR,BER_zf_dfe_u(uu,:),'bv:',...
             SNR,BER_mmse_dfe_u(uu,:),'rx:',SNR,BER_mmse_osic_u2(uu,:),'ks:')
    grid;
    legend('ZFE','MMSE-LE','ZF-DFE','MMSE-DFE','MMSE-OSIC','ZF-DFE-bef-decoupling','MMSE-DFE-bef-decoupling','MMSE-OSIC-bef-decoupling');
    % legend('ZFE','MMSE-LE','ZF-DFE','MMSE-DFE');%,'northeast');
    xlabel('SNR in dB');
    ylabel('Bit Error Rate');
    title('Performance comparison of different equalization techniques for N = 8 and nu = 2');
    % axis([5,25,10^(-6),10^0]);
end






% % figure(1);
% % set(gcf, 'Name','User 1'); % or set(figure(1), 'Name', 'User 1');
% % % figure(aaa);
% % % aaa = figure('Name', 'User1');
% % % semilogy(SNR,BER_zfe_1,'b*--',SNR,BER_mmse_le_1,'ms--',SNR,BER_zf_dfe_1,'ro:',SNR, BER_mmse_dfe_1,'y^--')
% % semilogy(SNR,BER_zfe(1,:),'b*--',SNR,BER_mmse_le(1,:),'ms--',SNR,BER_zf_dfe(1,:),'ro:',SNR,BER_zf_dfe_u1,'gv:',SNR, BER_mmse_dfe(1,:),'y^--')
% % grid;
% % legend('ZFE','MMSE-LE','ZF-DFE','ZF-DFE-U','MMSE-DFE');
% % % legend('ZFE','MMSE-LE','ZF-DFE','MMSE-DFE');%,'northeast');
% % xlabel('SNR in dB');
% % ylabel('Bit Error Rate');
% % title('Performance comparison of different equalization techniques for STBC in ISI multiple access channels');
% % % axis([5,25,10^(-6),10^0]);
% % 
% % figure(2);
% % set(gcf, 'Name','User 2');
% % % semilogy(SNR,BER_zfe_2,'b*--',SNR,BER_mmse_le_2,'ms--',SNR,BER_zf_dfe_2,'ro:',SNR, BER_mmse_dfe_2,'y^--')
% % semilogy(SNR,BER_zfe(2,:),'b*--',SNR,BER_mmse_le(2,:),'ms--',SNR,BER_zf_dfe(2,:),'ro:',SNR,BER_zf_dfe_u2,'gv:',SNR, BER_mmse_dfe(2,:),'y^--')
% % grid;
% % legend('ZFE','MMSE-LE','ZF-DFE','ZF-DFE-U','MMSE-DFE');
% % % legend('ZFE', 'MMSE-LE', 'ZF-DFE', 'MMSE-DFE');%, 'northeast');
% % xlabel('SNR in dB');
% % ylabel('Bit Error Rate');
% % title('Performance comparison of different equalization techniques for STBC in ISI multiple access channels');
% % 
% % figure(3);
% % set(gcf, 'Name','User 3');
% % % semilogy(SNR,BER_zfe_3,'b*--',SNR,BER_mmse_le_3,'ms--',SNR,BER_zf_dfe_3,'ro:',SNR, BER_mmse_dfe_3,'y^--')
% % semilogy(SNR,BER_zfe(3,:),'b*--',SNR,BER_mmse_le(3,:),'ms--',SNR,BER_zf_dfe(3,:),'ro:',SNR,BER_zf_dfe_u3,'gv:',SNR, BER_mmse_dfe(3,:),'y^--')
% % grid;
% % legend('ZFE','MMSE-LE','ZF-DFE','ZF-DFE-U','MMSE-DFE');
% % % legend('ZFE', 'MMSE-LE', 'ZF-DFE', 'MMSE-DFE');%, 'northeast');
% % xlabel('SNR in dB');
% % ylabel('Bit Error Rate');
% % title('Performance comparison of different equalization techniques for STBC in ISI multiple access channels');
