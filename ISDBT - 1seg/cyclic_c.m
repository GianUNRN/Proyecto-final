close all; clear; clc;

% n = 184  ; % Codeword length
% k = 102; % Message length
% 
% data = randi([0 1], 1, k); %RMSB X^3 + X


n=7;
k=4;
g = [1,1,0,1];
data = [1,0,0,0];

b = CyclycCode(data, n,k,g);

S_table = [ 1	0	0;
            0	1	0;
            0	0	1;
            1	1	0;
            0	1	1;
            1	1	1;
            1	0	1];

code = [b, data];

e = zeros(1,7);
e(randi([1 7])) = 1;

r = xor(code, e);

s = ParityCheck(r, g, n, k);

[~,idx] = ismember(s, S_table, 'rows');

correction = de2bi(2^(idx-1), n);

r_corr = xor(r, correction);


function parity_bits = CyclycCode(data, n, k, g)
    shift_data = [zeros(1, n-k), data];
    
    max_p_data = find(shift_data == 1);
    max_p_data = max_p_data(end);
    shift_data = shift_data(1:max_p_data);
    while max_p_data>=n-k+1
        
        divi = [zeros(1,max_p_data - (n-k)-1), g];
        
        r = xor(shift_data, divi);
        
        max_p_data = find(r == 1);
        if isempty(max_p_data)
            shift_data = zeros(1, n-k);
            max_p_data = 0;
        else
            max_p_data = max_p_data(end);
            
            shift_data = r(1:max_p_data);
        end
    end
    parity_bits = [double(shift_data), zeros(1, (n-k)-length(shift_data))];
end

function syndrom = ParityCheck(r, g, n, k)
    max_p_r = find(r == 1);
    max_p_r = max_p_r(end);
    r = r(1:max_p_r);
    s = r;
    while max_p_r>=n-k+1
    
        divi = [zeros(1,max_p_r - (n-k)-1), g];
    
        s = xor(s, divi);
    
        max_p_r = find(s == 1);
        if isempty(max_p_r)
            s = 0;
            max_p_r = 0;
        else
            max_p_r = max_p_r(end);
    
            s = s(1:max_p_r);
        end
    end
    syndrom = [double(s), zeros(1, n-k -length(s))];
end

