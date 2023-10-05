clear;
clc;
close all;
%step-1)givens
%a)Number of bits/SNR=10^6 bits
Number_of_bits = 1000;
%b)signal to noise ratio range=0 to 30dB with 2dB steps
Snr_range = 0:2:20; % SNR Range In DB

%step-2)generating random binary data vector
% Generate binary data vector
Data = randi([0 1], 1, Number_of_bits);

%ASK Modulation (OOK)
Modulated_data_ask = Data;

%FSK Modulation (PRK) 
Modulated_data_fsk = zeros(1, Number_of_bits);
for I=1:Number_of_bits 
    if Data(I) == 0 
    Modulated_data_fsk(I) =1; 
    else
    Modulated_data_fsk(I) = 3+1i; 
    end
end

% PSK Modulation (PRK)
Modulated_data_psk = pskmod(Data, 2);

% Create PSK Modulator Object 
Psk_mod = modem.pskmod('m', 2);

% Modulate Data Using PSK 
Modulated_data_psk_modem = Psk_mod.modulate(Data);

% BER Calculation 
Ber_ask = zeros(1, length(Snr_range));
Ber_fsk = zeros(1, length(Snr_range));
Ber_psk = zeros(1, length(Snr_range));
Ber_psk_modem = zeros(1, length(Snr_range));
Ber_pam_modem = zeros(1, length(Snr_range));

% Define Threshold Value 
Threshold = 0.5;

% Loop Over SNR Range 
for I = 1:length(Snr_range) 
    % Compute Noise Variance From SNR 
    Snr = Snr_range(I); 
    Noise_var = 1/(10^(Snr/10));

    %Step-3)Applying Noise To Bits
    % Add AWGN To Data
    Received_data_ask = Modulated_data_ask + sqrt(Noise_var)*randn(1, Number_of_bits);
    Received_data_fsk = Modulated_data_fsk + sqrt(Noise_var/2)*(randn(1, Number_of_bits)+1i*randn(1, Number_of_bits))*sqrt(2);
    Received_data_psk = Modulated_data_psk + sqrt(Noise_var/2)*randn(1,Number_of_bits);
    Received_data_psk_modem = awgn(Modulated_data_psk_modem, Snr, 'measured');

    %Step-4)Decidin Wheter Th Rx_Sequence Is 1 Or 0
    
    %Step-4)Demodulation
    % Demodulate Received Data
    Detected_bits_ask = (Received_data_ask > Threshold);
    Detected_bits_fsk = zeros(1, Number_of_bits);
    for J=1:Number_of_bits
        if real(Received_data_fsk(J)) > 2
        Detected_bits_fsk(J)=1;
        end
    end
    Detected_bits_psk = (Received_data_psk > 0);
    Detected_bits_psk_modem = pskdemod(Received_data_psk_modem, 2);

    %Step-5)Computing BER
    % Compute BER
    Ber_ask(I) = sum(xor(Data, Detected_bits_ask))/Number_of_bits;
    Ber_fsk(I) = sum(xor(Data, Detected_bits_fsk))/Number_of_bits;
    Ber_psk(I) = sum(xor(Data, Detected_bits_psk))/Number_of_bits;
    Ber_psk_modem(I) = sum(xor(Data, Detected_bits_psk_modem))/Number_of_bits;
end

% Plot Results Figure;
semilogy(Snr_range, Ber_ask, 'B-O'); hold on; 
semilogy(Snr_range, Ber_fsk, 'R-O');
semilogy(Snr_range, Ber_psk_modem, 'M-O'); 
grid on;
xlabel('SNR (DB)'); ylabel('Bit Error Rate (BER)'); 
legend('ASK (OOK)', 'FSK', 'PSK (BPSK)', 'PSK (Modem.Pskmod)', 'PAM (Modem.Pammod)');


