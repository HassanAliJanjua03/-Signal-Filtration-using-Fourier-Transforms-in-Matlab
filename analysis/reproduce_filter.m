% Reproducible frequency-domain analysis of the Signals.m experiment.
% Park-Miller + Box-Muller noise for matching the published numerical example.
% Original coursework script Signals.m remains unchanged.
clear; close all; clc;
fs = 1000; N = 1000; t = (0:N-1)/fs;
f1 = 50; f2 = 150; fc = 100;
desired = sin(2*pi*f1*t);
signal = desired + 0.5*sin(2*pi*f2*t);
state = 42; noise = zeros(1,N);
for k = 1:N
    state = mod(16807*state,2147483647); u1 = state/2147483647;
    state = mod(16807*state,2147483647); u2 = state/2147483647;
    noise(k) = 0.5*sqrt(-2*log(u1))*cos(2*pi*u2);
end
noisy = signal + noise;
frequencies = (-N/2:N/2-1)*fs/N;
X = fftshift(fft(noisy));
H = double(abs(frequencies)<=fc);
Y = X.*H;
filtered = real(ifft(ifftshift(Y)));
rmse = @(x) sqrt(mean((x-desired).^2));
snr_db = @(x) 10*log10(sum(desired.^2)/sum((x-desired).^2));
fprintf('Input RMSE: %.6f\nOutput RMSE: %.6f\n',rmse(noisy),rmse(filtered));
fprintf('Input SNR: %.6f dB\nOutput SNR: %.6f dB\n',snr_db(noisy),snr_db(filtered));
assert(max(abs(Y(abs(frequencies)>fc)))==0);
assert(max(abs(Y(abs(frequencies)<=fc)-X(abs(frequencies)<=fc)))==0);
figure('Name','Time-domain comparison');
subplot(3,1,1); plot(t,desired); xlim([0 .1]); ylim([-3 3]); grid on; title('50 Hz reference');
subplot(3,1,2); plot(t,noisy); xlim([0 .1]); ylim([-3 3]); grid on; title('Noisy input');
subplot(3,1,3); plot(t,filtered); xlim([0 .1]); ylim([-3 3]); grid on; title('Filtered output'); xlabel('Time (s)');
% One-sided amplitude spectra: DC and Nyquist are not doubled.
A = abs(fft(noisy))/N; B = abs(fft(filtered))/N;
A = A(1:N/2+1); B = B(1:N/2+1);
A(2:end-1)=2*A(2:end-1); B(2:end-1)=2*B(2:end-1);
figure('Name','Spectrum comparison');
plot((0:N/2)*fs/N,A,'Color',[.8 .38 .22]); hold on;
plot((0:N/2)*fs/N,B,'Color',[.41 .25 .65]); grid on;
xlabel('Frequency (Hz)'); ylabel('Amplitude'); legend('Noisy input','Filtered output');
figure('Name','Ideal low-pass mask'); stairs(frequencies,H,'LineWidth',2); grid on;
xlabel('Frequency (Hz)'); ylabel('Mask magnitude'); ylim([-.1 1.1]);
