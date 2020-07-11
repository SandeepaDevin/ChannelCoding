function [ Error_count, BER ] = Hamming_15_11( Bit_count, E_b, Noise_var )
%BPSK Summary of this function goes here
%   Detailed explanation goes here

H_Bit_count = Bit_count*15/11;
Tx_Bit_array = randi([0 1], Bit_count, 1);      % Bit array to be transmitted
Tx_Symbol_array = zeros(H_Bit_count, 1);        % Symbol array to be transmitted
Rx_Hamming_array = zeros(H_Bit_count, 1);       % Array to store bits from Rx symbols
Error_count = 0;

% Hamming Coder
Window_count = Bit_count/11;                    % Assumes length to be multiple of four
Eleven_bit_windows = reshape(Tx_Bit_array, 11, Window_count)';

P_Matrix = [   1 1 0 0;
               0 1 1 0;
               0 0 1 1;
               1 0 1 0;
               1 0 0 1;
               0 1 0 1;
               1 1 1 0;
               0 1 1 1;
               1 0 1 1;
               1 1 0 1;
               1 1 1 1; ];                       % Coeffcient Matrix
G_Matrix = [P_Matrix eye(11)];                   % Generator Matrix

Tx_Hamming_array = mod(Eleven_bit_windows * G_Matrix, 2);
Tx_Hamming_array = reshape(Tx_Hamming_array', 1, numel(Tx_Hamming_array))';

% Generating the signal array based on the bit array
for q = 1:H_Bit_count
    if Tx_Hamming_array(q) == 1
        Tx_Symbol_array(q) = sqrt(E_b);
    else 
        Tx_Symbol_array(q) = -sqrt(E_b);
    end
end

% Generating the received signal by adding signal and noise arrays
Noise_array = wgn(H_Bit_count, 1, 10*log10(Noise_var), 'complex');
Rx_Symbol_array = Tx_Symbol_array + Noise_array;

% Making the decision based on decision boundary
for q = 1:H_Bit_count
    if abs(Rx_Symbol_array(q)-1) <= abs(Rx_Symbol_array(q)+1)
        Rx_Hamming_array(q) = 1;
    else
        Rx_Hamming_array(q) = 0;
    end
end

% Hamming Decoder
Window_count = numel(Rx_Hamming_array)/15; 
Fifteen_bit_windows = reshape(Rx_Hamming_array, 15, Window_count)';

H_Matrix = [eye(4) P_Matrix'];                  % Parity Check Matrix
S_Table = syndtable(H_Matrix);                  % Syndrome Table
Syndrome = mod(Fifteen_bit_windows*H_Matrix',2);
S_Table_index = Syndrome*[8; 4; 2; 1];

E_Matrix = zeros(Window_count, 15);              % Error Matrix
for q = 1:Window_count
    E_Matrix(q,:) = S_Table(S_Table_index(q)+1,:);
end

Rx_Hamming_array = mod(Fifteen_bit_windows + E_Matrix, 2); % Performing error correction
Rx_Hamming_array(:,1:4) = [];                            % Removing parity bits
Rx_Bit_array = reshape(Rx_Hamming_array', Bit_count, 1);

% Error Counting
for q = 1:Bit_count
    if Rx_Bit_array(q)~=Tx_Bit_array(q)
        Error_count = Error_count+1;        
    end
end

BER = Error_count/Bit_count;

end

