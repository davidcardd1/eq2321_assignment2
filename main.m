
load('assignment2.mat');

Fs = 8000;
y = speech8;
% The Uniform Scalar Quantization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = -6:0.01:6;
% xmax = 4;
% n_bits = 2;
% m=0;
% idx = sq_enc(x,n_bits,xmax,m);
% outq = sq_dec(idx,n_bits,xmax,m);
% 
% figure(1);
% subplot(1,2,1);
% plot(x,outq);
% grid("on");
% legend("m=0");
% 
% m=1.5;
% idx = sq_enc(x,n_bits,xmax,m);
% outq = sq_dec(idx,n_bits,xmax,m);
% subplot(1,2,2);
% plot(x,outq);
% legend("m=1.5");
% grid("on");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parametric Coding of Speech%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[E, V, A, P] = analysis(y, 256, 8, 10);

% figure(2);
h_lin = histogram(E);hold on;
m_lin = mean(h_lin.BinEdges);
xmax_lin = max(h_lin.BinEdges) - m_lin;
xline(m_lin,'-r',{'m'});
xline(m_lin+xmax_lin,'-g',{'xmax'});
xline(m_lin-xmax_lin,'-g',{'xmax'});
n_bits_lin = 6;
idx = sq_enc(E,n_bits_lin, xmax_lin, m_lin);
outq = sq_dec(idx,n_bits_lin, xmax_lin, m_lin);

% figure(3);
% plot(E); 
% hold on; 
% plot(outq);
% legend('Energy', 'Quantized linear-Energy')

[x_lin] = synthesis(outq, V, A, P, 8);

% log
E_log = log(E);
% figure(4)
h_log = histogram(E_log);

m_log = mean(h_log.BinEdges);
xmax_log = max(h_log.BinEdges) - m_log;
xline(m_log,'-r',{'m'});
xline(m_log+xmax_log,'-g',{'xmax'});
xline(m_log-xmax_log,'-g',{'xmax'});
n_bits_log = 3;

idx = sq_enc(E_log, n_bits_log, xmax_log, m_log);
outq_log = sq_dec(idx, n_bits_log, xmax_log, m_log);

% figure(5)
% plot(E); 
% hold on; 
% plot(exp(outq_log));
% legend('Energy', 'Quantized log-Energy')

[x_log] = synthesis(exp(outq_log), V, A, P, 8); 
% soundsc(x_lin, Fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quantizing the Pitch and Voiced/Unvoiced Decision%%%%%%%%%%%%%%%%

P_log = log(P);
% figure(6)
h_pitch_log = histogram(P_log);

m_pitch_log = mean(h_pitch_log.BinEdges);
xmax_pitch_log = max(h_pitch_log.BinEdges) - m_pitch_log;
xline(m_pitch_log,'-r',{'m'});
xline(m_pitch_log+xmax_pitch_log,'-g',{'xmax'});
xline(m_pitch_log-xmax_pitch_log,'-g',{'xmax'});
n_bits_pitch_log = 5;

idx = sq_enc(P_log, n_bits_pitch_log, xmax_pitch_log, m_pitch_log);
outq_plog = sq_dec(idx, n_bits_pitch_log, xmax_pitch_log, m_pitch_log);

% figure(7)
plot(P); hold on; plot(exp(outq_plog));
legend('Pitch', 'Quantized log-Pitch')

[x_pitch_log] = synthesis(E, V, A, exp(outq_plog), 8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quantizing the LP parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
codeA = encodefilter(A, lsfCB1, lsfCB2);
Aq = decodefilter(codeA, lsfCB1, lsfCB2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimizing the Bit Allocation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xq] = synthesis(exp(outq_log), V, Aq, exp(outq_plog), 8);

SNR = 10*log10(var(y(1:length(xq)))/cov(y(1:length(xq)) - xq));

bpsp = (n_bits_log + n_bits_pitch_log + 1 + 20)/8; % bit per sample
bps = (n_bits_log + n_bits_pitch_log + 1 + 20)*Fs/8; % bit per second

bps
bpsp





