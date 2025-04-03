close all; clear; clc;


%%
addpath('funciones\')
load("c_data.mat")

%%
bit_intrl = Coder.Bit_Intrlv(coded_data);

symb_intrl = Coder.Symb_Intrlv(bit_intrl).';

ak = Coder.Bi2QAM(symb_intrl);
 
%%



wk = Coder.wk;

k_continual_pilot = [0 48 54 87 141 156 192 201 255 279 282 333 432 450 483 525 531 618 636 714 759 765 780 804 873 888 918 939 942 969 984 1050 1101 1107 1110 1137 1140 1146 1206 1269 1323 1377 1491 1683 1704];
continual_p = wk(k_continual_pilot+1);

%%
TPS_K = [34 50 209 346 413 569 595 688 790 901 1073 1219 1262 1286 1469 1594 1687];

s_0 = 1-2*wk(TPS_K + 1);





%%


ofdm_frame = zeros(68,1705);
Tps = Coder.TpsBits(21,1,0, 3/4,0, 2, 64, 1/4);


TPS_Matrix = zeros(68, length(TPS_K));
TPS_Matrix(1,:) = s_0.';
TPS_Matrix(2:end,:) = TPS_Matrix(2:end,:) + Tps;


for j= 2:68
    TPS_Matrix(j,:) = (TPS_Matrix(j,:)*(-2) + 1).*TPS_Matrix(j-1,:);
end

ak_sin_uso = ak;

for i = 1:68
    ofdm_frame(i,k_continual_pilot+1) = 4/3 *(1 - 2*continual_p);
    
    k_scatter_pilot = Coder.Sct_Pilots(i-1, 1704);
    scatter_p = wk(k_scatter_pilot+1);
    ofdm_frame(i, k_scatter_pilot+1) = 4/3 *(1 - 2*scatter_p);


    ofdm_frame(i,TPS_K+1) = TPS_Matrix(i,:);

    data_cells = find(ofdm_frame(i,:) == 0);
    num_data_cells = length(data_cells);

    if(length(ak_sin_uso) > num_data_cells)
        ofdm_frame(i, data_cells) = ak_sin_uso(1:num_data_cells);
        ak_sin_uso = ak_sin_uso(num_data_cells+1:end);
    end

    
end

clear j i;


%%

ofdm_frame = [zeros(68,171), ofdm_frame,zeros(68,172)];

ofdm_frame_time = ifft(ofdm_frame, 2048, 2);

[spec, ~] = pwelch(ofdm_frame_time(1,:), ones(1,2048), 0, 2048);
semilogy( spec)


%%

guard_p = 1/4;
signal = [ofdm_frame_time(:, end - 2048*guard_p +1:end), ofdm_frame_time];

signalrow = reshape(signal.', 1,[]);

%%

r = fft(ofdm_frame_time,2048,2);
r = r(:, 172:172+1704);

%%

aux = Decoder.Exctract_TPS(r);
dec = comm.BCHDecoder(67, 53, 'X^14 + X^9 + X^8 + X^6 + X^5 + X^4 + X^2 + X + 1');

tpb = dec(aux);
%%

idx_sync = Decoder.Find_Sync(tpb);

ofdm_sym = Decoder.Exctract_Pilots(r);
qam_sym = Decoder.Serial_Data(ofdm_sym);
%%
bits = Decoder.QAM2Bi(qam_sym);


