close all; clear; clc;
tf = 5;
Fs = 2000;          % Freq muestreo
Fc = 10e3;           % portadora
Fm = 1e3;            % f carrier
beta = 5;            % indice modulación
f_doppler = 0;

t = linspace(0, tf, tf*Fs );     
m = cos(2*pi*Fm*t);  

y_fm = cos(2*pi*Fc*t + beta*sin(2*pi*Fm*t));

[psi, freq, range] = caf(y_fm, y_fm .* exp(1i*2*pi*f_doppler*t), 100, Fs*tf, Fs);


imagesc(freq, range, abs(psi));
colormap('hot')
xlabel('Frequency Shift (Hz)');
ylabel('Range (m)');
function [psi, freq_axis, range_axis] = caf(xe, xr, Np, N, fs)
    c = 3e8;
    xe_b = reshape(xe, [Np, N/Np]);
    xr_b = reshape(xr, [Np, N/Np]);
    
    D = ifft(fft(xe_b,2*Np-1,1) .* conj(fft(xr_b,2*Np-1,1)));
    D = fftshift(D, 1);
    psi = fft(D,[],2);
    
    psi = fftshift(psi,2);
    
    freq_axis = linspace(-fs*(2*Np - 1)/2, fs*(2*Np - 1)/2, size(psi, 2));
    
    range_resolution = (1 / fs) * c / 2;
    range_axis = linspace(-size(psi, 1)/2, size(psi, 1)/2, size(psi, 2))* range_resolution;
end