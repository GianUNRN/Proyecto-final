close all; clear; clc;
N = 2470; 

Np = 19;
c = 3*10^8;
n = 0:N-1;
f_doppler = 3e4; 
fs = 10e6;
delay_samples = 3; 



xr = randn(1,N);

xe_doppler = xr.* exp(1j * 2 * pi * f_doppler * n / fs);

clutt = [zeros(1, 10) xr(1:end-10)] +  [zeros(1, 3) xr(1:end-3)] + [zeros(1, 6) xr(1:end-6)];
xe = [zeros(1, delay_samples) xe_doppler(1:end-delay_samples)] + clutt;


xe_f = eca_clutter_filter(xr,xe, 1000, N).';
figure(1)
[CAF, f, r] = caf(xe_f, xr, Np, N, fs);
imagesc(f, r, abs(CAF));
colormap('hot')

figure(2)
[CAF_2, f_2, r_2] = caf(xe, xr, Np, N, fs);
imagesc(f_2, r_2, abs(CAF_2));
colormap('hot')


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
    
    freq_axis = linspace(-fs*(2*Np - 1)/2, fs*(2*Np - 1)/2, size(psi, 2));
    
    range_resolution = (1 / fs) * c / 2;
    range_axis = linspace(-size(psi, 1)/2, size(psi, 1)/2, size(psi, 2))* range_resolution;
end