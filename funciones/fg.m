classdef fg
   methods(Static)
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
        
        function Tps = TpsBits(numbits, frame, alpha, Rate_Hp, Rate_Lp, Modo, Orden_Const, t_guard)
            
                    sync_tps = [[0,0,1,1,0,1,0,1,1,1,1,0,1,1,1,0]; 
                                [1,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1]];
                
                    length_tps = [[0,1,0,1,1,1]; % Cell Identification information is not transmitted (23 TPS bits in use)
                                  [0,1,1,1,1,1]]; %Cell Identification information is transmitted (31 TPS bits in use)
                
                
                    frame_tps = [[0,0];     %frame 1
                                 [0,1];     %frame 2
                                 [1,0];     %frame 3
                                 [1,1]];    %frame 4
                
                
                    const  = [[0,0];    %QPSK
                              [0,1];    %16 QAM
                              [1,0]];   % 64 QAM
                
                
                    hierarchy_tps = [[0,0,0];   %sin jerarquia
                                     [0,0,1];   %alpha = 1
                                     [0,1,0];   %alpha = 2
                                     [0,1,1]];  %alpha = 4
                    
                    ratesHP_tps =  [[0,0,0];   %1/2
                                    [0,0,1];   %2/3
                                    [0,1,0];   %3/4
                                    [0,1,1];   %5/6
                                    [1,0,0]];  %7/8
                
                    ratesLP_tps =  [[0,0,0];   %1/2
                                    [0,0,1];   %2/3
                                    [0,1,0];   %3/4
                                    [0,1,1];   %5/6
                                    [1,0,0]];  %7/8
                
                    guard_intr = [[0,0];    %1/32
                                  [0,1];    %1/16
                                  [1,0];    %1/8
                                  [1,1]];   %1/4
                    
                    transmision_tps = [[0,0];    %2k
                                       [0,1]];   %8k
                
                    Tps = NaN;
                    
                    
                    %sync
                    if frame == 1 || frame == 3
                        aux = sync_tps(1,:);
                    elseif frame == 2 || frame == 4
                        aux = sync_tps(2,:);
                    else
                        return
                    end
                    
                    %length
                    if numbits == 21
                        aux = [aux, length_tps(1,:)];
                    elseif numbits == 31
                        aux = [aux, length_tps(2,:)];
                    else
                        return
                    end
                
                    %frame number
                    if 1<=frame && frame<=4
                        aux = [aux, frame_tps(frame,:)];
                    else
                        return
                    end
                
                    %constelation
                    if Orden_Const == 4
                        aux = [aux, const(1,:)];
                    elseif Orden_Const == 16
                        aux = [aux, const(2,:)];
                    elseif Orden_Const == 64
                        aux = [aux, const(3,:)];
                    else
                        return
                    end
                
                    
                    
                    %hierarchy
                    if alpha == 0
                        aux = [aux, hierarchy_tps(1,:)];
                    elseif alpha == 1
                        aux = [aux, hierarchy_tps(2,:)];
                    elseif alpha ==2
                        aux = [aux, hierarchy_tps(3,:)];
                    elseif alpha ==4
                        aux = [aux, hierarchy_tps(4,:)];
                    else
                        return
                    end
                
                
                    %Data Rate Hp
                    if Rate_Hp == 1/2
                        aux = [aux, ratesHP_tps(1,:)];
                    elseif Rate_Hp == 2/3
                        aux = [aux, ratesHP_tps(2,:)];
                    elseif Rate_Hp == 3/4
                        aux = [aux, ratesHP_tps(3,:)];
                    elseif Rate_Hp == 5/6
                        aux = [aux, ratesHP_tps(4,:)];
                    elseif Rate_Hp == 7/8
                        aux = [aux, ratesHP_tps(5,:)];
                    else
                        return
                    end
                    
                    %Data Rate Lp
                    if alpha == 0
                        aux = [aux, zeros(1,3)];
                    else
                        if Rate_Lp == 1/2
                            aux = [aux, ratesLP_tps(1,:)];
                        elseif Rate_Lp == 2/3
                            aux = [aux, ratesLP_tps(2,:)];
                        elseif Rate_Lp == 3/4
                            aux = [aux, ratesLP_tps(3,:)];
                        elseif Rate_Lp == 5/6
                            aux = [aux, ratesLP_tps(4,:)];
                        elseif Rate_Lp == 7/8
                            aux = [aux, ratesLP_tps(5,:)];
                        else
                            return
                        end
                    end
                    %Guard interval
                    if t_guard == 1/36
                        aux = [aux, guard_intr(1,:)];
                    elseif t_guard == 1/16
                        aux = [aux, guard_intr(2,:)];
                    elseif t_guard == 1/8
                        aux = [aux, guard_intr(3,:)];
                    elseif t_guard == 1/4
                        aux = [aux, guard_intr(4,:)];
                    else 
                        return
                    end
                
                    
                
                    %Transmision mode
                    if Modo == 2
                        aux = [aux, transmision_tps(1,:)];
                    elseif Modo == 8
                        aux = [aux, transmision_tps(2,:)];
                    else 
                        return
                    end
                
                    %Cell ID
                    aux = [aux, zeros(1,8)];
                
                    %zeros
                    aux = [aux, zeros(1,6)];
                    
                    enc = comm.BCHEncoder(67, 53, 'X^14 + X^9 + X^8 + X^6 + X^5 + X^4 + X^2 + X + 1');
                    Tps = enc(aux.');
                
                
                
                end
        function Gen = PRBSGen(Kmax)
            Gen = comm.PNSequence('Polynomial', [11 2 0], ...  % x^11 + x^2 + 1
                                  'InitialConditions', ones(1,11), ...  % Seed
                                  'SamplesPerFrame', Kmax+1);
        end
        function k_idx = Sct_Pilots(l, Kmax)
            in = 3*mod(l,4);
            k_idx = in:12:Kmax;
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


        function tps = Exctract_TPS(frame)

            TPS_K = [34 50 209 346 413 569 595 688 790 901 1073 1219 1262 1286 1469 1594 1687];
            TPS_symb = frame(:, TPS_K+1);
            tps = -1*ones(67,1);
            for j=2:68

                if mean(abs(TPS_symb(j,:) - TPS_symb(j-1,:))) < 1e-3
                    tps(j-1) = 0;
                else 
                    tps(j-1) = 1;

                end
            end

    
        end
    

   end
end
