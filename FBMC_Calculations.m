
%snrdB = 12;              % SNR in dB
SNR=1:15;
for snrdB=1:length(SNR)
s = rng(211);            % Set RNG state for repeatability
numFFT = 512;           % Number of FFT points
numGuards = 212;         % Guard bands on both sides
K = 4;                   % Overlapping symbols, one of 2, 3, or 4
numSymbols = 100;        % Simulation length in symbols

% Prototype filter
switch K
    case 2
        HkOneSided = sqrt(2)/2;
    case 3
        HkOneSided = [0.911438 0.411438];
    case 4
        HkOneSided = [0.971960 sqrt(2)/2 0.235147];
    otherwise
        return
end
% Build symmetric filter
Hk = [fliplr(HkOneSided) 1 HkOneSided];

% QAM symbol mapper
bitsPerSubCarrier = 2;
qamMapper = comm.RectangularQAMModulator(...
    'ModulationOrder', 2^bitsPerSubCarrier, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% Transmit-end processing
%   Initialize arrays
L = numFFT-2*numGuards;  % Number of complex symbols per OFDM symbol
KF = K*numFFT;
KL = K*L;
dataSubCar = zeros(L, 1);
dataSubCarUp = zeros(KL, 1);

sumFBMCSpec = zeros(KF*2, 1);
sumOFDMSpec = zeros(numFFT*2, 1);

numBits = bitsPerSubCarrier*L/2;    % account for oversampling by 2
inpData = zeros(numBits, numSymbols);
rxBits = zeros(numBits, numSymbols);
txSigAll = complex(zeros(KF, numSymbols));
symBuf = complex(zeros(2*KF, 1));

for S = [1:1:30]
sefbmc = (bitsPerSubCarrier * S )/ (S + K - 0.5);
end




% Loop over symbols
for symIdx = 1:numSymbols

    % Generate mapped symbol data
    inpData(:, symIdx) = randi([0 1], numBits, 1);
    modData = step (qamMapper, inpData(:, symIdx));

    % OQAM Modulator: alternate real and imaginary parts
    if rem(symIdx,2)==1     % Odd symbols
        dataSubCar(1:2:L) = real(modData);
        dataSubCar(2:2:L) = 1i*imag(modData);
    else                    % Even symbols
        dataSubCar(1:2:L) = 1i*imag(modData);
        dataSubCar(2:2:L) = real(modData);
    end

    % Upsample by K, pad with guards, and filter with the prototype filter
    dataSubCarUp(1:K:end) = dataSubCar;
    dataBitsUpPad = [zeros(numGuards*K,1); dataSubCarUp; zeros(numGuards*K,1)];
    X1 = filter(Hk, 1, dataBitsUpPad);
    % Remove 1/2 filter length delay
    X = [X1(K:end); zeros(K-1,1)];

    % Compute IFFT of length KF for the transmitted symbol
    txSymb = fftshift(ifft(X));

    % Transmitted signal is a sum of the delayed real, imag symbols
    symBuf = [symBuf(numFFT/2+1:end); complex(zeros(numFFT/2,1))];
    symBuf(KF+(1:KF)) = symBuf(KF+(1:KF)) + txSymb;

    % Compute power spectral density (PSD)
    currSym = complex(symBuf(1:KF));
    [specFBMC, fFBMC] = periodogram(currSym, hann(KF, 'periodic'), KF*2, 1);
    sumFBMCSpec = sumFBMCSpec + specFBMC;

    % Store transmitted signals for all symbols
    txSigAll(:,symIdx) = currSym;
end


%% OFDM Modulation with Corresponding Parameters
%
% For comparison, we review the existing OFDM modulation technique, using
% the full occupied band, however, without a cyclic prefix.

for symIdx = 1:numSymbols
    
    inpData2 = randi([0 1], bitsPerSubCarrier*L, 1);
    modData = step (qamMapper, inpData2);
        
    symOFDM = [zeros(numGuards,1); modData; zeros(numGuards,1)];
    ifftOut = sqrt(numFFT).*ifft(ifftshift(symOFDM));

    [specOFDM,fOFDM] = periodogram(ifftOut, rectwin(length(ifftOut)), ...
        numFFT*2, 1, 'centered'); 
    sumOFDMSpec = sumOFDMSpec + specOFDM;
end




% QAM demodulator
qamDemod = comm.RectangularQAMDemodulator(...
    'ModulationOrder', 2^bitsPerSubCarrier, ...
    'BitOutput', true, ...
    'NormalizationMethod', 'Average power');
BER = comm.ErrorRate;

% Process symbol-wise
for symIdx = 1:numSymbols
    rxSig = txSigAll(:, symIdx);

    % Add WGN
    rxNsig = awgn(rxSig, snrdB, 'measured');

    % Perform FFT
    rxf = fft(fftshift(rxNsig));

    % Matched filtering with prototype filter
    rxfmf = filter(Hk, 1, rxf);
    % Remove K-1 delay elements
    rxfmf = [rxfmf(K:end); zeros(K-1,1)];
    % Remove guards
    rxfmfg = rxfmf(numGuards*K+1:end-numGuards*K);

    % OQAM post-processing
    %  Downsample by 2K, extract real and imaginary parts
    if rem(symIdx, 2)
        % Imaginary part is K samples after real one
        r1 = real(rxfmfg(1:2*K:end));
        r2 = imag(rxfmfg(K+1:2*K:end));
        rcomb = complex(r1, r2);
    else
        % Real part is K samples after imaginary one
        r1 = imag(rxfmfg(1:2*K:end));
        r2 = real(rxfmfg(K+1:2*K:end));
        rcomb = complex(r2, r1);
    end

    % Demapper: Perform hard decision
    rxBits(:, symIdx) = step(qamDemod, rcomb);
end

% Measure BER with appropriate delay
BER.ReceiveDelay = 2*KL;
ber = step (BER, inpData(:), rxBits(:));

% Display Bit error
disp(['FBMC Reception for K = ' num2str(K) ', BER = ' num2str(ber(1)) ...
    ' at SNR = ' num2str(snrdB) ' dB'])

% Restore RNG state
rng(s);
a(snrdB)=ber(1);
disp (a(snrdB));
end
% Plot power spectral density (PSD) over all subcarriers
sumOFDMSpec = sumOFDMSpec/mean(sumOFDMSpec(1+2*numGuards:end-2*numGuards));
figure; 
plot(fOFDM,10*log10(sumOFDMSpec)); 
hold on
grid on
axis([-0.5 0.5 -180 10]);
xlabel('Normalized frequency'); 
ylabel('PSD (dBW/Hz)')
title(['Power Spectral Density - FBMC vs OFDM']);

% Plot power spectral density
sumFBMCSpec = sumFBMCSpec/mean(sumFBMCSpec(1+K+2*numGuards*K:end-2*numGuards*K-K));
plot(fFBMC-0.5,10*log10(sumFBMCSpec));
grid on
axis([-0.5 0.5 -180 10]);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['Power Spectral Density - FBMC vs OFDM']);

figure
semilogy (SNR,a, '--xg');
axis ([1 15 0 0.005 ]);
xlabel('Signal to Noise Ratio in dB');
ylabel('Bit Error Rate');
title ('SNR vs BER - FBMC');
grid on;
hold off;

% Compute peak-to-average-power ratio (PAPR)
PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprFBMC] = step (PAPR, txSymb);
disp(['Peak-to-Average-Power-Ratio (PAPR) for FBMC = ' num2str(paprFBMC) ' dB']);

% Compute peak-to-average-power ratio (PAPR)
PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprOFDM] = step(PAPR2,ifftOut);
disp(['Peak-to-Average-Power-Ratio (PAPR) for OFDM = ' num2str(paprOFDM) ' dB']);