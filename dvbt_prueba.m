close all; clear; clc;
addpath('funciones\')
cargar_data = 0;
Coder = Coder(1/2,2,64);
Decoder = Decoder(2);
Decoder = Decoder.Set_Mod(64);
Decoder = Decoder.Set_Rate(1/2);
%%
if cargar_data
    load("c_data.mat"); %#ok<*NOPRT> 
    
else
    
    rng(23); %#ok<*NOPRT> 
    sym = randi([0 255],1, 188*800);
    coded_data = Coder.CD(sym);
end


%%

[bit_intrl, left_bits] = Coder.Bit_Intrlv(coded_data);


[symb_intrl, left_symbs] = Coder.Symb_Intrlv(bit_intrl);

ak = Coder.Bi2QAM(symb_intrl.');
 

[ofdm_frame, ak_sin_uso] = Coder.Frame(ak,21,1,0, 3/4, 0, 1/4);




%%

ofdm_frame_pad = [zeros(68,Coder.pad), ofdm_frame,zeros(68,Coder.pad+1)];


ofdm_frame_time = ifft(ofdm_frame_pad, Coder.F_max, 2);


[spec, ~] = pwelch(ofdm_frame_time(1,:), ones(1,Coder.F_max), 0, Coder.F_max);
semilogy( spec)


% %%
% 
% guard_p = 1/4;
% signal = [ofdm_frame_time(:, end - 2048*guard_p +1:end), ofdm_frame_time];
% signal2 =  [ofdm_frame_time2(:, end - 2048*guard_p +1:end), ofdm_frame_time2];
% 
% signalrow = reshape(signal.', 1,[]);
% signalrow2 = reshape(signal.', 1,[]);
% 
% signal_r = [signalrow, signalrow2];
% 
% [cor, lags] = xcorr( signal_r(1:600),signal_r(1:20000));

%%

r = fft(ofdm_frame_time,Coder.F_max,2);
r = r(:, Coder.pad+1:Coder.pad +1 + Coder.K_max);

%%

aux = Decoder.Exctract_TPS(r);
dec = comm.BCHDecoder(67, 53, 'X^14 + X^9 + X^8 + X^6 + X^5 + X^4 + X^2 + X + 1');

tpb = dec(aux);


%%

idx_sync = Decoder.Find_Sync(tpb);

ofdm_sym = Decoder.Exctract_Pilots(r); 


qam_sym = Decoder.Serial_Data(ofdm_sym);

%%
r_intrl_symb = Decoder.QAM2Bi(qam_sym);

r_bits = Decoder.Symb_Deintrlv(r_intrl_symb.');

r_data = Decoder.Bits_Deintrlv(r_bits);

%%

r_dec = Decoder.Viterbi_Dec(r_data);


r_deintrl = Decoder.Conv_Deintrlv(r_dec);
%%

[dec_data, data_left] = Decoder.RS_Deco(r_deintrl);


display(isequal(dec_data, sym(1:length(dec_data))))


%%

aux = 1:2000;
aux2 = 4000:6000;
[intrlvData, state] =convintrlv(aux, Coder.nrows, Coder.slope);

[intrlvData2, state2] =convintrlv(aux2, Coder.nrows, Coder.slope,state);
p2= convdeintrlv(intrlvData2,Coder.nrows, Coder.slope);




