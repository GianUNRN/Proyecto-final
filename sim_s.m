close all; clear; clc;
N = 2470; 

Np = 19;
c = 3*10^8;
n = 0:N-1;
f_doppler = 7e4; 
fs = 10e6;
delay_samples = 6; 



xe = randn(1,N);

xr_doppler = xe.* exp(1j * 2 * pi * f_doppler * n / fs);


xr = [zeros(1, delay_samples) xr_doppler(1:end-delay_samples)];

xe_b = reshape(xe, [Np, N/Np]);
xr_b = reshape(xr, [Np, N/Np]);

D = ifft(fft(xe_b,2*Np-1,1) .* conj(fft(xr_b,2*Np-1,1)));
D = fftshift(D, 1);
psi = fft(D,[],2);

psi = fftshift(psi,2);




figure(1)




freq_axis = linspace(-fs*(2*Np - 1)/2, fs*(2*Np - 1)/2, size(psi, 2));

range_resolution = (1 / fs) * c / 2;
range_axis = linspace(-size(psi, 1)/2, size(psi, 1)/2, size(psi, 2))* range_resolution;

imagesc(freq_axis, range_axis, abs(psi));
colormap('hot')
xlabel('Frequency Shift (Hz)');
ylabel('Range (m)');
