
classdef PR
   methods(Static)
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
    
   end
end
