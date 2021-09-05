%this function returns the bits after decoding the received signals.
%received_signal is a Mt by 1 matrix of complex numbers. Mt is the number
%of transmitted antenna or it can be any number. For each qpsk symbols it
%generates 2 bits and then it makes a vector of all the bits of the
%received symbols. So the returned vector will have a length of 1 by 2Mt

function x = qpsk_decoded_bits(received_signal)

%this is a 2 by length(received_signal) matrix. 1 is added with each
%element and then each of them is divided by 2 to make the -1 to 0.
bits = [sign(real(received_signal))+1 sign(imag(received_signal))+1]./2;

x = [];

for k = 1:length(received_signal)
    x = [x bits(k,:)];
end


% %if received_signal is a vector of 1 by Mt elements
% 
% bits = [sign(real(received_signal))+1 sign(imag(received_signal))+1]./2;
% 
% x = [];
% 
% len = length(received_signal);
% for k = 1:len
%     x = [x bits(k) bits(k+len)];
% end
