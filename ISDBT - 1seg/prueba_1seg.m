close all; clear; clc;

load("data_mpeg.mat")
data_mpeg = data_mpeg(1:131600);
%%
Rate = 2/3;
Const = 4;

Coder = Coder(Rate,8,Const);
Decoder = Decoder(8);
Decoder = Decoder.Set_Mod(Const);
Decoder = Decoder.Set_Rate(Rate);
%%

cod_data = Coder.RS(data_mpeg);
sym = Coder.RandMpeg(cod_data);
%%
sym_pad= [sym;zeros(Coder.D,1)];
intrlv_bytes = Coder.Conv_Intrl(sym_pad);
intrlvbits = de2bi(intrlv_bytes);
intrlvbits = reshape(intrlvbits', 1, []);

x = Coder.Conv_Code(intrlvbits)';

%%
bit_intrl = Coder.Bit_Intrlv(x);

ak = Coder.Bi2QAM(bit_intrl);
%%
div = floor(length(ak)/384)*384;

ak_data_seg = reshape(ak(1:div), 384,[]).';

permuted_ak(:,Coder.Freq_intrl_per_m3+1) = ak_data_seg;
AC = randi([0, 1], 1,203);
capa_A = [[1, 1, 0, 1]; 
          [1, 1, 0, 1]]; 

capa_B = [[7, 7, 7, 15]; 
          [7, 7, 7, 15]]; 
capa_C = [[7, 7, 7, 15]; 
          [7, 7, 7, 15]];

Frame = Coder.Frame(permuted_ak,AC,0,1, 1, 1, 1, capa_A, capa_B, capa_C);
%%
Ofdm_frame = [zeros(204,2592), Frame, zeros(204,2593)];

Ofdm_frame_time = ifft(ifftshift(Ofdm_frame), 2^13,2);
Ofdm_frame_time = [Ofdm_frame_time(:,end-256 + 1 :end), Ofdm_frame_time];
%%
tx_signal = reshape(Ofdm_frame_time.', [], 1);
[cor, lags] = xcorr(tx_signal);
plot(lags, abs(cor))

