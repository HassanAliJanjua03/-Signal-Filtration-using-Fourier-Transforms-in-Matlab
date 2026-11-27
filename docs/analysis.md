# Numerical analysis

The original experiment in Signals.m uses a 50 Hz sine wave, a 150 Hz sine wave at amplitude 0.5, and Gaussian noise with standard deviation 0.5. The variable called signal in that script contains both sinusoidal components. For this reproducible comparison the desired reference is the 50 Hz sine wave alone.

## Reproduction parameters

- Sample rate: 1000 Hz, sample count: 1000, duration: 1 second.
- FFT bin spacing: 1 Hz; centered bins: -500 through 499 Hz.
- Ideal mask: 1 for absolute frequency at most 100 Hz, 0 otherwise.
- Noise sequence: Park–Miller state initialized to 42, multiplier 16807 and modulus 2147483647. Each sample consumes two uniforms via the Box–Muller cosine formula. This differs from the original unseeded randn call.
- DFT and inverse DFT were evaluated independently to generate the published plots and metrics. The MATLAB FFT reproduction has not been run in MATLAB.

## Numerical results

| Quantity | Value |
| --- | ---: |
| Input RMSE | 0.589268 |
| Output RMSE | 0.234397 |
| Input SNR (dB) | 1.583450 |
| Output SNR (dB) | 9.590645 |
| SNR gain (dB) | 8.007195 |
| Input 50 Hz amplitude | 1.014523 |
| Input 150 Hz amplitude | 0.463432 |
| Output 50 Hz amplitude | 1.014523 |
| Output 150 Hz amplitude | 0.000000 |

RMSE is computed relative to the desired 50 Hz sinusoid, and SNR uses the ratio of reference energy to error energy. The one-sided spectrum doubles non-DC and non-Nyquist magnitudes after dividing by N. Time-domain graphics show the first 100 ms of the full record. The original project report is retained unchanged; these figures are derived from the MATLAB source, not extracted from the report.

## Limitations

The mask retains in-band noise and removes out-of-band signal components regardless of their origin. An ideal rectangular frequency response is not a causal real-time filter. FFT-domain filtering assumes a periodic extension of the finite record; boundary effects and ringing can occur. The SNR improvement depends on the noise realization and is not an experimental hardware measurement.

[Back to project](../README.md)
