function C= WIENER_HOPF_MMSE(vsampled ,param, OP , chdata, txffe, Noise_XC,ivs) 
cmx=param.RxFFE_cmx;
cpx=param.RxFFE_cpx;
% do this early on so we can reuse the old code
% to be replaced with MMSE function from healey...
if param.N_bg ~=0 % must be floating taps
    cpx=param.N_bmax; % N_f in spreadsheet
end
num_taps=cmx+cpx+1;
ndfe=param.ndfe;
spui=param.samples_per_ui;
%% Start of WIENER-HOPF MMSE EQ code
                R_n = zeros(num_taps,num_taps);
                R_xt = zeros(num_taps,num_taps);

                if OP.Do_Colored_Noise
                    % Form Noise Autocorrelation matrix
                    Noise_XC = reshape (Noise_XC,1,[]);
                    len = length(Noise_XC);
                    if len < num_taps, Noise_XC = [Noise_XC,zeros(1,num_taps-len)]; end
                    Noise_XC = Noise_XC(1:num_taps);
                	R_n = toeplitz (Noise_XC,Noise_XC);
                end
                %% Calculate Cross Talk Correlation matrix at T intervals.
                if OP.Do_XT_Noise
                    % Calculate variance of Tx signal based on +/-1 outer limits
                    Tx_sigma = sqrt(  (param.levels^2-1)/(3*(param.levels-1)^2) );
                    for jj = 2:length(chdata)
                        if isequal(chdata(jj).type, 'FEXT') || isequal(chdata(jj).type, 'NEXT')
                            if isequal(chdata(jj).type, 'FEXT')
                                % len = length(chdata(jj).ctle_imp_response);
                                % ch_imp = ifft( fft(reshape(chdata(jj).ctle_imp_response,1,[])) .* fft(reshape(upsample(txffe, param.samples_per_ui),1,[]),len) );
                                ch_imp= FFE( txffe , cmx,param.samples_per_ui, chdata(jj).ctle_imp_response).';
                                ch_imp = filter(ones(param.samples_per_ui, 1), 1, ch_imp);
                            elseif isequal(chdata(jj).type, 'NEXT')
                                ch_imp = chdata(jj).ctle_imp_response;
                                ch_imp = filter(ones(param.samples_per_ui, 1), 1, ch_imp);
                            end
                            norms = zeros(1,spui);
                            for ii = 1:spui
                                norms(ii) = norm(ch_imp(ii:spui:end));
                            end
                            % Pick out sampling phase with largest noise contribution
                            [~,cursor] = max(norms);
                            sub_sample_ch = ch_imp(cursor:spui:end);
                            xc = xcorr(sub_sample_ch,num_taps)* Tx_sigma^2;
                            xc = xc(num_taps+1:end);
                            xc = xc(1:num_taps);
                            R = toeplitz (xc,xc);
                            R_xt = R_xt + R;
                        end
                    end
                end
                %% Noise + Cross Talk contribution to R matrix
                R_n_xc = zeros(num_taps+ndfe);
                R_n_xc(1:num_taps,1:num_taps) = R_n + R_xt ;
                %% For least means squares, we want to solve
                %
                %  ro =  |Ryy  Ryx| * w
                %        |Rxy  Rxx|
                % see Cioffi chapter 3, 3.7.3

                himp = vsampled;
                RefTap = cmx+1;
                Signal_Variance = mean ( [-1, -1/3, 1/3, 1].^2 );
                Ryy.r = [1:num_taps];
                Ryy.c = [1:num_taps];
                Ryx.r  = 1:num_taps;
                Ryx.c  = num_taps + (1:ndfe);
                Rxy.r  = num_taps + (1:ndfe);
                Rxy.c  = 1:num_taps;
                Rxx.r  = num_taps+(1:ndfe);
                Rxx.c  = num_taps+(1:ndfe);
                himp = reshape (himp,1,[]);
                %% ro is simply the channel response reversed in time
                himp_lr = fliplr(himp); % make sure himp has enough pre/post cursors, e.g. numtaps+ndfe
                himp_lr = [zeros(1,num_taps+ndfe), himp_lr, zeros(1,num_taps+ndfe)];
                [~,pk] = max(himp_lr);
                r_indx = (1:num_taps) - RefTap;
                ro = [himp_lr(pk+r_indx),  zeros(1,ndfe)].';
                ro = ro*Signal_Variance;
                %% Setup up the covariance matrix
                R = zeros(num_taps+ndfe);
                % Form Ryy
                % Note: important to use whole impulse response
                % not just the part that spans the FFE.
                [ryy,lags] = xcorr(himp,himp, num_taps-1);
                R(Ryy.r,Ryy.c) = toeplitz (ryy(num_taps:end));

                % Form Rxx
                R(Rxx.r,Rxx.c) = diag(ones(1,ndfe));

                % Form Ryx columns
                Ryx_indxs = (1:num_taps)-1;
                for jj = 0:ndfe-1
                    %             R(Ryx.r,Ryx.c(jj+1) ) = himp_lr(pk+1 + Ryx_indxs -jj ).';
                    R(Ryx.r,Ryx.c(jj+1) ) = himp_lr(pk-cmx-1 + Ryx_indxs -jj ).';
                end
                % Form Rxy rows
                R(Rxy.r,Rxy.c ) = R(Ryx.r,Ryx.c )';

                % add in Signal Variance
                R = R*Signal_Variance;
                Rtmp = R;
                % Add in Xt and colored noise terms
                R = R + R_n_xc;

                % SNR = 25 dB
                SNR = OP.FFE_SNR;
                Noise_var =   max(vsampled)^2*Signal_Variance/10^(SNR/10);
                R_noise = diag(ones(1,num_taps))*Noise_var;
                if ~OP.Do_Colored_Noise &&  ~OP.Do_XT_Noise
                    R(Ryy.r,Ryy.c) = R(Ryy.r,Ryy.c) + R_noise;
                end


                %% Solve for equalizer weights
                w = inv(R)*ro;
                C = w;
                %% Deal with 1st post Cursor DFE weight saturation
                %  ro = Rw by moving "saturated" weights over to the LHS
                DFE_h1_indx = num_taps+1;
                Indx_full = 1:length(C);
                ws = C;
                if ndfe>0 && abs(C(DFE_h1_indx)) > param.bmax(1)
                    rtmp = reshape (ro,[],1);
                    Rtmp = R;
              	    % Move saturated DFE weights over to left hand side of equation
                    ws = zeros (size(C));
                    ws (DFE_h1_indx) = sign(C(DFE_h1_indx))*param.bmax(1);
             	    rtmp = rtmp - Rtmp*ws;

             	    % and remove the corresponding column from R
              	    Rtmp(:,DFE_h1_indx) = [];
                    Indx_full (DFE_h1_indx) = [];
                    % now Rtmp isn't square so have to use the R'R trick
                    % Probably a little dicey "theoretically" because
                    % w = inv(R)*ro is already the mmse solution
                    % now we at doing a R'R operation, but hey
                    % seems to work
                    % Alternative, since R is now over specified, more rows than
                    % columns, one could try removing one of the DFE rows from the
                    % Rxy  Rxx portion of the R matrix.

              	    w_partial = inv(Rtmp'*Rtmp)*(Rtmp'*rtmp);
              	    ws (Indx_full,:) = w_partial;
                    C = ws;
                end
                % From Cioffi, Chapter 3
                var_ffe_dfe = Signal_Variance-ws.'*ro;
                SNR_Eq = 10*log10(Signal_Variance/var_ffe_dfe-1);

                %% Scale FFE gain to target output voltage
                Target_ouput = vsampled(ivs)*10^(param.current_ffegain/20);
                C = C*Target_ouput;
                C = C(1:num_taps);
            %% End MMSE dfe code