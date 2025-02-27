close all; clear; clc;

format long
N = 2500; 

Np = 50;
c = 3*10^8;
n = 0:N-1;
f_doppler = 2.5e4; 
fs = 10e6;
delay_samples = 3; 

noise = rand(1, N);


xr = randn(1,N);

xe_doppler = xr.* exp(1j * 2 * pi * f_doppler * n / fs);

clutt = xr;
for i = 1:200
    clutt = clutt + [zeros(1, i) xr(1:end-i)];
end

xe = ([zeros(1, delay_samples) xe_doppler(1:end-delay_samples)]*0.5 + clutt);


figure(1)
[CAF, f, r] = caf(xe, xr, Np, N, fs);
imagesc(f, r, abs(CAF));
colormap('hot')

title('Señal original')

%-----------------------------------------------------------------------------%
figure(2)

xe_f = eca_clutter_filter(xr,xe, 1000, N);

[CAF_2, f_2, r_2] = caf(xe_f, xr, Np, N, fs);
imagesc(f_2, r_2, abs(CAF_2));
xticks(linspace(min(f_2), max(f_2), length(f_2)))
yticks(linspace(min(r_2), max(r_2), length(r_2)))
colormap('hot')
colorbar()
title('Filtro LS')

%-----------------------------------------------------------------------------%
%{
tic()
xe_f_p = block_lattice_filter(xr,xe, 1000);
toc()
figure(3)
[CAF_3, f_3, r_3] = caf(xe_f_p, xr, Np, N, fs);
imagesc(f_3, r_3, abs(CAF_3));
colormap('hot')
title('filtro lattice')


function e = block_lattice_filter(xr, xe, M)

    N = length(xr);
    b = zeros(M + 1, N); 
    f = zeros(M + 1, N);
    e = zeros(M + 1, N); 
    h = zeros(1, M + 1); 
    k = zeros(1, M); 
    
    b(1, :) = xr; 
    f(1, :) = xr; 
    e(1, :) = xe; 
  
    for i = 1:M
        num = 2 * sum(b(i, 1:end-1) .* f(i, 2:end));
        den = sum(abs(f(i, 2:end)).^2 + abs(b(i, 1:end-1)).^2);
        k(i) = num / den;
        
    
        b(i+1, 2:end) = b(i, 1:end-1) - k(i) * conj(f(i, 2:end));
        f(i+1, 1:end-1) = f(i, 2:end) - conj(k(i)) * b(i, 1:end-1);

        h(i) = sum(e(i, :) .* conj(b(i, :))) / sum(abs(b(i, :)).^2);
        e(i+1, :) = e(i, :) - h(i) * b(i, :);
    end

    e = e(end, :).';

end
%}
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

