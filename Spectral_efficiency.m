grid on;
bitsPerSubCarrier=4;
K=4;
filterLen=43;
NCP=43;
numFFT=1024;
S=[1:2:30]
%fbmc
z=S*bitsPerSubCarrier*10;
i=(S+K-0.5)*10;
seFBMC=z./i;
semilogy(S,seFBMC,'--*b')
%ufmc
hold on
a=bitsPerSubCarrier*numFFT*S;
b=(numFFT+filterLen-1)*S;
seUFMC=a./b;
semilogy(S,seUFMC,'--Or')
%ofdm
hold on
q=bitsPerSubCarrier*numFFT*S;
r=(numFFT+NCP)*S;
seOFDM=q./r;
semilogy(S,seOFDM,'--*y')


legend('FBMC','UFMC','OFDM','Location','SouthEast')
    title('Spectral Efficiency-FBMC vs UFMC vs OFDM')
    xlabel('Duration of Burst[ms]')
    ylabel('Spectral Efficiency [b/s/Hz] for m=4')
    grid on
    hold off
return ;