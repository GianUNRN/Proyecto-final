close all; clear; clc;

%%
coded_data = CD(567);

%%
bit_intrl = Bit_Intrlv(coded_data);
symb_intrl = Symb_Intrlv(bit_intrl);

%%
qpsk_symb = bi2de(symb_intrl','left-msb');
a_k = qammod(qpsk_symb, 2^size(symb_intrl,1), 'gray' );
%%
function x = CD(m)
    x = zeros(2176*m,1);

    n= 204;
    k= 188;
    p_gen = 'D8 + D4 + D3 + D2 + 1';
    N = 188;

    % Interleaver Convolucional
    nrows = 12; %cant de shift registers
    slope = 1;  %diferencia de delays entre ramas

    % Codigo Convolucional con perforaciones
        
    code_gen = [171 133]; 
    len = 7; 
    
    trellis = poly2trellis(len,code_gen);
    
    punc_pat = [1;0;1;1;1;0]; % R = 3/4 ==> X1 Y1 Y2 X3
    for i = 0:m-1
        
        sym = randi([0 255],1, N);
        
        % Reed Solomon
        
        
        msg = gf(sym, 8, p_gen);
        y = rsenc(msg, n, k);
        
        
        
        intrlvData =convintrlv(y.x, nrows, slope);
        
        intrlvbits = de2bi(intrlvData)';
        intrlvbits = intrlvbits(:);
        
        
        coded_data = convenc(intrlvbits, trellis, punc_pat);

        x(2176*i+1: 2176*i+2176) =  coded_data;
    end
    
end


function bit_intrl = Bit_Intrlv(data)
    bit_intrl = NaN;
    if mod(data, 756) ~= 0
        return;
    end
    demux_data = reshape(data, 6, []); 

    H_bit = zeros(6, 126);
    % Funciones de permutacion
    H_bit(1, :) = (0:125) + 1;
    H_bit(2, :) = mod((0:125) + 63, 126) + 1;
    H_bit(3, :) = mod((0:125) + 105, 126)+ 1;
    H_bit(4, :) = mod((0:125) + 42, 126)+ 1;
    H_bit(5, :) = mod((0:125) + 21, 126)+ 1;
    H_bit(6, :) = mod((0:125) + 84, 126)+ 1;
    
    m = size(demux_data,2)/126;
    
    idx = zeros(size(H_bit,1), m*size(H_bit,2));
    
    for i = 0:m-1
        idx(:, 126*i+1: 126*i+126) =  H_bit + 126*i;
    end
    bit_intrl = zeros(size(idx));
    bit_intrl(1,:) = demux_data(1,idx(1,:));
    bit_intrl(2,:) = demux_data(2,idx(2,:));
    bit_intrl(3,:) = demux_data(3,idx(3,:));
    bit_intrl(4,:) = demux_data(4,idx(4,:));
    bit_intrl(5,:) = demux_data(5,idx(5,:));
    bit_intrl(6,:) = demux_data(6,idx(6,:));
end

function symb_intrl = Symb_Intrlv(data)
    N_max = 1512; % 2k mode
    symb_intrl = NaN;
    if mod(size(data,2), N_max) ~= 0
        return;
    end
   
    M_max = 2048;
    Nr = log2(M_max);
    Rp = zeros(M_max-1, Nr - 1);
    
    Rp(1,:) = 0;
    Rp(2,:) = 0;
    Rp(3,:) = [1 zeros(1, Nr - 2)];
    
    for j = 3:M_max-1
        Rp(j+1,1:end-1) =  Rp(j,2:end);
        Rp(j+1,end) = xor(Rp(j,1),Rp(j,4));
    end
    
    R_permutation = [4, 3, 9, 6, 2, 8, 1, 5, 7, 0] + 1;
    
    R = Rp(:, R_permutation);
    
    q = 0;
    H = zeros(1, N_max);
    for i = 1: M_max
        aux = R(i, :) * (2.^(0:Nr-2)).';
    
        H(q+1) = mod(i-1, 2)* 2^(Nr-1) + aux;
    
        if(H(q+1) < N_max)
            q = q + 1;
        end
    end
    H = H + 1;
    m = size(data,2)/N_max;
    idx = zeros(1, m*length(H));
    
    for i = 0:m-1
        idx(:, N_max*i+1: N_max*i+N_max) =  H + N_max*i;
    end

    symb_intrl = 3*zeros(size(data));
    
    for k = 0:2:size(data,2)-1
        symb_intrl(:,idx(k+1)) = data(:, k+1); 
        symb_intrl(:,k+2) = data(:, idx(k+2)); 
    end
    


end