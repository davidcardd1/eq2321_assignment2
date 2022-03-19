
load('assignment2.mat');

Fs = 8000;

%% The Uniform Scalar Quantization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = -6:0.01:6;
xmax = 4;
n_bits = 2;
m=0;
idx = sq_enc(x,n_bits,xmax,m);
outq = sq_dec(idx,n_bits,xmax,m);

figure(1);
subplot(1,2,1);
plot(x,outq);
grid on;
title("m=0");

m=1.5;
idx = sq_enc(x,n_bits,xmax,m);
outq = sq_dec(idx,n_bits,xmax,m);
subplot(1,2,2);
plot(x,outq);
title("m=1.5");
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parametric Coding of Speech

% Quantizing the Gain%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = speech8;
[E, V, A, P] = analysis(x, 256, 8, 10);

% task 1
figure();
histL = histogram(E);hold on;
meanL = mean(histL.BinEdges);
xmaxL = max(histL.BinEdges) - meanL;
xline(meanL,'-r',{'m =', meanL});
xline(meanL+xmaxL,'-g',{'m+xmax = ', meanL+xmaxL});
xline(meanL-xmaxL,'-g',{'m-xmax =' , meanL-xmaxL});

%task 2
n_bitsL = 6; % try/error with this parameter
idx = sq_enc(E,n_bitsL, xmaxL, meanL);
outq = sq_dec(idx,n_bitsL, xmaxL, meanL);
soundsc(synthesis(outq, V, A, P, 8), fs);

% % figure(3);
% % plot(E); 
% % hold on; 
% % plot(outq);
% % legend('Energy', 'Quantized linear-Energy')
% 
% [x_lin] = synthesis(outq, V, A, P, 8);

%task 3
figure();
histPLog = histogram(log(E));hold on;
meanLog = mean(histPLog.BinEdges);
xmaxLog = max(histPLog.BinEdges) - meanLog;
xline(meanLog,'-r',{'m =', meanLog});
xline(meanLog+xmaxLog,'-g',{'m+xmax = ', meanLog+xmaxLog});
xline(meanLog-xmaxLog,'-g',{'m-xmax =' , meanLog-xmaxLog});

%task 4
n_bitsLog = 5; %try/error
idx = sq_enc(log(E), n_bitsLog, xmaxLog, meanLog);
outqLog = sq_dec(idx, n_bitsLog, xmaxLog, meanLog);
soundsc(synthesis(exp(outqLog), V, A, P, 8), fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quantizing the Pitch and Voiced/Unvoiced Decision%%%%%%%%%%%%%%%%

%pitch
histPLog = histogram(log(P));hold on;
meanPLog = mean(histPLog.BinEdges);
xmaxPLog = max(histPLog.BinEdges) - meanPLog;
n_bitsPLog = 6; %try/error
idx = sq_enc(log(P), n_bitsPLog, xmaxPLog, meanPLog);
outqPLog = sq_dec(idx, n_bitsPLog, xmaxPLog, meanPLog);
soundsc(synthesis(E, V, A, exp(outqPLog), 8), fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quantizing the LP parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
codeA = encodefilter(A, lsfCB1, lsfCB2);
Aq = decodefilter(codeA, lsfCB1, lsfCB2);

soundsc(synthesis(E, V, Aq, P, 8), fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimizing the Bit Allocation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%task 1
xhat = synthesis(exp(outqLog), V, Aq, exp(outqPLog), 8);
xhat2 = synthesis(E, V, Aq, P, 8);
soundsc(xhat, fs)
SNR = 10*log10(var(x(1:length(xhat2)))/var(x(1:length(xhat2)) - xhat2))

bpsp = (n_bitsLog + n_bitsPLog + 1 + 20)/8; % bit per sample
bps = (n_bitsLog + n_bitsPLog + 1 + 20)*Fs/8; % bit per second

bps
bpsp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Speech Waveform Quantization

% Uniform Scalar Quantization of Speech %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Adaptive Open-Loop DPCM


