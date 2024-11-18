close all; clear; clc;
%----------------------------------------------------------------------------------------%
%% Version bruta y literal de CAF
N = 40; 

xe = randn(1, N)*exp(1j*2*pi*2); 
xr = xe; 

psi = zeros(2*N+1, 2*N+1);

for m = 1:2*N+1
    for k = 1:2*N+1    
        for n = 0:N-1
            if(n-(m-1-N)>= 0 && n-(m-1-N) <= N-1)         
                psi(m, k) = psi(m, k) + xe(n+1) * conj(xr(n-(m-1-N)+1)) * exp(-1j * 2 * pi * (k-1-N) * n / N);
            end
        end        
    end    
end
figure(1)
R = -N:N;
colormap('hot');
imagesc(R,R,abs(psi));
colorbar

%----------------------------------------------------------------------------------------%
%% Primer intento alg Batches
N = 40; 

xe = randn(1, N); 
xr = xe; 

figure(2)
Np = 4;
xe_b = reshape(xe, [Np, N/Np]);
xr_b = reshape(xr, [Np, N/Np]);

D = zeros(2*Np - 1, N/Np);
tic
for j = 1:N/Np 
    D(:,j) = xcorr(xe_b(:,j),xr_b(:,j));
end


psi_b = fft(D,[],2);
psi_b = fftshift(psi_b,2);
toc
colormap('hot');
imagesc(abs(psi_b));
colorbar

%----------------------------------------------------------------------------------------%
%% Segundo intento alg Batches
tic
D_aux_f = ifft(fft(xe_b,2*Np-1,1) .* conj(fft(xr_b,2*Np-1,1)));
D_aux_f = fftshift(D_aux_f, 1);
psi_b = fft(D_aux_f,[],2);
psi_b = fftshift(psi_b,2);
toc

%más óptimo, aprobecha mas la complejidad de fft / ifft

figure(3)
colormap('hot');
imagesc(abs(psi_b));
colorbar