close all; clear; clc;
%%
addpath('funciones\')
cargar_data = 1;

if cargar_data
    load("c_data.mat")%#ok
else
    rng(23); %#ok
    sym = randi([0 255],1, 188*3);
    coded_data = Coder.CD(sym);
end


%%

bit_intrl = Coder.Bit_Intrlv(coded_data);


symb_intrl = Coder.Symb_Intrlv(bit_intrl).';

ak = Coder.Bi2QAM(symb_intrl);
 

Tps = Coder.TpsBits(21,1,0, 3/4,0, 2, 64, 1/4);

[ofdm_frame, ak_sin_uso] = Coder.Frame(ak, Tps);

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

% idx_sync = Decoder.Find_Sync(tpb);

ofdm_sym = Decoder.Exctract_Pilots(r, 1704);
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
