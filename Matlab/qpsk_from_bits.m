    function X = qpsk_from_bits(D)

    % Angle [pi/4 3*pi/4 -3*pi/4 -pi/4] corresponds to 
    % Gray code vector [00 10 11 01], respectively.
   
    full_len = length (D);
    
    table = exp(j*[-3/4*pi 3/4*pi 1/4*pi -1/4*pi]);  % generates QPSK symbols on four different phase
    %table = exp(j*[-pi pi -2*pi 2*pi]);
    

    table = table([0 1 3 2]+1); % Gray code mapping pattern for QPSK symbols i.e. 00, 10, 11, 01
    
    bit_table = reshape(D,2,full_len/2);
    
    % maps transmitted bits into QPSK symbols. every qpsk symbol is divided
    % by 2 because there are four transmit antennas and the power should be
    % normalized so that every complex qpsk symbol has 1/4 energy and they
    % sum upto unit energy at the receiver
    X = (table([2 1]*bit_table+1));  
    