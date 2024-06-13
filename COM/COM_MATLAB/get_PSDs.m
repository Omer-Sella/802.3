function result = get_PSDs(result,h,cursor_i, txffe,G_DC,G_DC2,param,chdata,OP)
if 1 % force indent for doc
    num_ui=param.num_ui_RXFF_noise;
    M=param.samples_per_ui;
    L=param.levels;
    sigma_X2=(L^2-1)/(3*(L-1)^2);
    f_b=param.fb;
    T_b=1/f_b;
    delta_f = f_b/num_ui; % Units are Hz.
    fvec = (0:num_ui*M/2)*delta_f; % Single-sided frequency axis.
    result.fvec=fvec;
    SNR_TX=param.SNR_TX;
    eta_0=param.eta_0; %V^2/GHz
end
if OP.WO_TXFFE % to speed up loop find sn onlu first time ctle is case
    %% compute S_rn healey_3dj_01_2401 slide 5
    S_RN_of_f=S_RN(fvec,G_DC,G_DC2,param);
    rxn_psd=[real(S_RN_of_f(1)), S_RN_of_f(2:end-1), real(S_RN_of_f(end)), conj(S_RN_of_f(end-1:-1:2))]; % Convert single-sided frequency response to conjugate-symmetric
    rxn_psd=rxn_psd/1e9;% Units are V^2/Hz.
    rxn_rms = sqrt(sum(rxn_psd)* delta_f);
    S_rn = sum(reshape(rxn_psd, num_ui, M).');
    S_rn=S_rn(1:num_ui/2+1);
    S_rn= [real( S_rn(1)),  S_rn(2:end-1), real( S_rn(end)), conj( S_rn(end-1:-1:2))];
    %rxn_acf_samp = ifft(S_rn)*f_b; % Autocorrelation function. Note that the first term is rxn_rms^2.
    if 0 % for debug
        figure
        set(gcf, 'tag', 'COM');movegui(gcf,'northeast');
        plot(fvec(1:num_ui)/f_b,10*log10((S_rn)*1000/100) ...
            ,'disp','Srn')
        xlim([0 0.5])
        ylim([-200 -140])
        set(gcf,'defaulttextinterpreter','none')
        xlabel('Normalized Frequency')
        ylabel('PSD dBm/Hz')
        hold on
        grid on
        title('PSD')
    end
    result.S_rn=S_rn;
    result.S_rn_rms=rxn_rms;

else % find noise for item that set have tx ffe for each loop
    %% from healey_3dj_01_2401 slide 6
    % Crosstalk power spectral density
    result.S_xn=0;
    if length(chdata)~=1;
        for xchan=2:length(chdata)
            pulse_ctle=filter(ones(1,M),1,chdata(xchan).ctle_imp_response(:).');
            pulse_ctle=[ pulse_ctle(1:floor(length(pulse_ctle)/M)*M) ];
            % enable less UI for computation speed improvement
            %%
            if num_ui*M > length(pulse_ctle)
                pulse_ctle= [ pulse_ctle zeros(1,num_ui*M-length(pulse_ctle)) ];
            else
                pulse_ctle=pulse_ctle(1:num_ui*M);
            end
            cmx=find(txffe==max(txffe))-1;
            for i1=1:M
                if ~strcmp(chdata(xchan).type,'NEXT')
                    hk(xchan).k=FFE( txffe , cmx,M, pulse_ctle )'; % need to speed up here
                else
                    hk(xchan).k=pulse_ctle;
                end
            end
            for i1=1:M
                hxn(i1)=norm(hk(xchan).k(i1:M:length(hk(xchan).k)) );
            end
            iphase(xchan)=find(hxn==max(hxn));
            hk(xchan).hrn= hk(xchan).k(iphase(xchan):M:length(hk(xchan).k));
            result.hk(xchan).hrn= hk(xchan).hrn;
            hk(xchan).S_xn=sigma_X2*(abs(fft(hk(xchan).hrn))).^2/param.fb;
            result.S_xn=hk(xchan).S_xn+result.S_xn;
        end
        result.xn_rms = sqrt(sum(result.S_xn)* delta_f);
    else
        result.xn_rms=0;
        iphase=0;
    end
    %% from healey_3dj_01_2401 slide 7
    % Transmitter noise power spectral density
    htn=filter(ones(1,M),1,chdata(1).ctle_imp_response); % ctle_imp_response does not have TxFFE included
    htn=htn(mod(cursor_i,M)+1:end-mod(cursor_i,M)); % align to sample point
    htn=reshape(htn,1,[]); % make row vectors
    htn=[ htn(1:floor(length(htn)/M)*M) ];
    htn= [htn zeros(1,num_ui*M-length(htn)) ]; 
    htn=htn(1:M:end);% resample
    if num_ui>length(htn)
        hext=[htn zeros(1,num_ui-length(htn))];
    else
        hext=htn(1:num_ui);
    end
    result.S_tn=sigma_X2*10^(-SNR_TX/10)*(abs(fft(hext))).^2/param.fb; % this corresponds to +/- pi
    result.tn_rms = sqrt(sum(result.S_tn)* delta_f); 
    %% from healey_3dj_01_2401 slide 8
    % Power spectral density of noise due to jitter
    %% Eq. 93A-28 %%
    sampling_offset = mod(cursor_i, M);
    %ensure we can take early sample
    if sampling_offset<=1
        sampling_offset=sampling_offset+M;
    end
    if (OP.LIMIT_JITTER_CONTRIB_TO_DFE_SPAN)
        cursors_early_sample = h(cursor_i-1+M*(-1:param.ndfe));
        cursors_late_sample = h(cursor_i+1+M*(-1:param.ndfe));
    else
        cursors_early_sample = h(sampling_offset-1:M:end);
        cursors_late_sample = h(sampling_offset+1:M:end);
    end
    % ensure lengths are equal
    cursors_early_sample = cursors_early_sample(1:length(cursors_late_sample));
    h_J = (cursors_late_sample-cursors_early_sample)/2*M;
    h_J=reshape(h_J,1,[]); % make row vectors
    if num_ui>length(h_J)
        h_J=[h_J zeros(1,num_ui-length(h_J))];
    else
        h_J=h_J(1:num_ui);
    end
    result.iphase=iphase;
    result.S_jn=sigma_X2*(param.A_DD^2+param.sigma_RJ^2)*(abs(fft(h_J))).^2/param.fb; % this corresponds to +/- pi
    result.jn_rms = sqrt(sum(result.S_jn)* delta_f);
    result.S_n=result.S_rn+ result.S_tn+ result.S_xn+ result.S_jn;
end