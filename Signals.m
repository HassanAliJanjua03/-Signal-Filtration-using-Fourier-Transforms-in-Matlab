% Define parameters
fs = 1000; % Sampling frequency
t = 0:1/fs:1-1/fs; % Time vector

f1 = 50; % Signal frequency
f2 = 150; % Noise frequency

% Create the signal with noise
signal = sin(2*pi*f1*t) + 0.5*sin(2*pi*f2*t); 
noise = 0.5 * randn(size(t)); 
noisy_signal = signal + noise; 

% Plot the original and noisy signals
figure; 
subplot(2, 1, 1); 
plot(t, signal, 'b', 'LineWidth', 1.5); 
title('Original Signal'); 
subplot(2, 1, 2); 
plot(t, noisy_signal, 'r', 'LineWidth', 1.5); 
title('Noisy Signal'); 

% Apply Fourier transform to analyze the frequency components
N = length(noisy_signal); 
frequencies = (-N/2:N/2-1)*(fs/N); 
fft_result = fftshift(fft(noisy_signal, N)); 
fft_magnitude = abs(fft_result); 

% Plot the magnitude spectrum
figure; 
plot(frequencies, fft_magnitude, 'b', 'LineWidth', 1.5); 
title('Frequency Spectrum of Noisy Signal'); 
xlabel('Frequency (Hz)'); 
ylabel('Magnitude'); 

% Design a low-pass filter in the frequency domain
cutoff_frequency = 100; % Adjust as needed
low_pass_filter = zeros(1, N); 
low_pass_filter(abs(frequencies) <= cutoff_frequency) = 1; 

% Apply the filter in the frequency domain
filtered_fft = fft_result .* low_pass_filter; 

% Plot the magnitude spectrum of the filtered signal
figure; 
plot(frequencies, abs(filtered_fft), 'r', 'LineWidth', 1.5); 
title('Frequency Spectrum of Filtered Signal'); 
xlabel('Frequency (Hz)'); 
ylabel('Magnitude'); 

% Apply inverse Fourier transform to obtain the filtered signal in the time domain
filtered_signal = ifft(ifftshift(filtered_fft)); 

% Plot the original, noisy, and filtered signals
figure; 
subplot(3, 1, 1); 
plot(t, signal, 'b', 'LineWidth', 1.5); 
title('Original Signal'); 
subplot(3, 1, 2); 
plot(t, noisy_signal, 'r', 'LineWidth', 1.5); 
title('Noisy Signal'); 
subplot(3, 1, 3); 
plot(t, real(filtered_signal), 'g', 'LineWidth', 1.5); 
title('Filtered Signal'); 
legend('Original Signal', 'Noisy Signal', 'Filtered Signal');
