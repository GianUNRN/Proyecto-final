close all; clear; clc;

N = 100; 

xe = randn(1, N); 
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

colormap('hot');
imagesc(abs(psi));
colorbar
%}