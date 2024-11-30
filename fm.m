close all; clear; clc;

fc = 100000;
kf = 5000; 

[audio, fs] = audioread('prueba.wav');
audio = mean(audio, 2); 





audio = audio / max(abs(audio));


int_audio = cumsum(audio) / fs; 


t = (0:length(audio)-1) / fs; 
fm_signal = cos(2 * pi * fc * t + 2 * pi * kf * int_audio.');
fm_signal = fm_signal(1:2500);


N = 2500;

Np = 50;
c = 3*10^8;
n = 0:N-1;
f_doppler = 200; 

delay_samples = 3; 

noise = rand(1, N);


xr = fm_signal + noise;

xe_doppler = xr.* exp(1j * 2 * pi * f_doppler * n / fs);

clutt = xr;
for i = 1:200
    clutt = clutt + [zeros(1, i) xr(1:end-i)];
end

xe = ([zeros(1, delay_samples) xe_doppler(1:end-delay_samples)]*0.5 + clutt + noise);


figure(1)
[CAF, f, r] = caf(xe, xr, Np, N, fs);
imagesc(f, r, abs(CAF));
colormap('hot')

figure(2)
tic()
xe_f = eca_clutter_filter(xr,xe, 1000, N);
toc()
[CAF_2, f_2, r_2] = caf(xe_f, xr, Np, N, fs);
imagesc(f_2, r_2, abs(CAF_2));
yticks(linspace(min(r_2), max(r_2), 20))
colormap('hot')
title('Filtro LS')



function y = eca_clutter_filter(xr, xe, K, N)
    
    xr = xr(:);
    xe = xe(:);

    R = zeros(N, K);
    for k = 1:K
        delay = k - 1;
        R(delay+1:end, k) = xr(1:N-delay); 
    end

    C = (R' * R) \ (R' * xe);

    clutter = R * C;

    y = xe - clutter;

    

end

function [psi, freq_axis, range_axis] = caf(xe, xr, Np, N, fs)

    c = 3e8;
    xe_b = reshape(xe, [Np, N/Np]);
    xr_b = reshape(xr, [Np, N/Np]);
    
    D = ifft(fft(xe_b,2*Np-1,1) .* conj(fft(xr_b,2*Np-1,1)));
    D = fftshift(D, 1);
    psi = fft(D,[],2);
    
    psi = fftshift(psi,2);
    
    freq_axis = linspace(-fs/(2*Np), fs/(2*Np), size(psi, 2));
    
    range_resolution = (1 / fs) * c;
    range_axis = linspace(-size(psi, 1)/2, size(psi, 1)/2, size(psi, 2))* range_resolution;
end
