close all;  clear; clc;
N = 188;
N0 = 1;

sym = randi([0 255],1, N);

%% Reed Solomon
n= 204;
k= 188;
p_gen = 'D8 + D4 + D3 + D2 + 1';

msg = gf(sym, 8, p_gen);
y = rsenc(msg, n, k);

%% Interleaver Convolucional
nrows = 12; %cant de shift registers
slope = 1;  %diferencia de delays entre ramas

intrlvData =convintrlv(y.x, nrows, slope);

intrlvbits = de2bi(intrlvData)';
intrlvbits = intrlvbits(:);


%% Codigo Convolucional con perforaciones

code_gen = [171 133]; 
len = 7; 

trellis = poly2trellis(len,code_gen);

punc_pat = [1;0;1;1;1;0]; % R = 3/4 ==> X1 Y1 Y2 X3

coded_data = convenc(intrlvbits, trellis, punc_pat);

%% Bit Interleaver 

%Demux
n_zeros = ceil(length(coded_data) / 756) * 756 - length(coded_data);

padd_data = [coded_data; zeros(n_zeros,1)];
demux_data = reshape(padd_data, 6, []); 

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

bit_intrl(1,:) = demux_data(1,idx(1,:));
bit_intrl(2,:) = demux_data(2,idx(2,:));
bit_intrl(3,:) = demux_data(3,idx(3,:));
bit_intrl(4,:) = demux_data(4,idx(4,:));
bit_intrl(5,:) = demux_data(5,idx(5,:));
bit_intrl(6,:) = demux_data(6,idx(6,:));
 
%% Symbol Interleaver 

N_max = 1512; % 2k mode
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

