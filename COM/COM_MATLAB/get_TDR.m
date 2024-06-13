function TDR_results = get_TDR(sdd, OP, param,ZT,np)
% sdd is differential s-parameters structure (2 port assumed)
% input parameter structure for s parameters sdd--> sdd.Impedance, sdd.Frequencies, sdd.Parameters, sdd.NumPorts
% TDR_results.delay             pre t=0 delay for TDR... help with time domain responce quaility
% TDR_results.tdr               the TDR responce (ohms vs  TDR_results.t
% TDR_results.t                 starting at t=0
% TDR_results.tx_filter         transmitter filter vs TDR_results.f
% TDR_results.Rx_filter         receiver filter vs TDR_results.f
% TDR_results.f                 frequency for filter and s parameters
% TDR_results.ptdr_RL           reflection waveform from the pulse
% TDR_results.WC_ptdr_samples_t worst case time sample of the reflection pulse
% TDR_results.WC_ptdr_samples   worst case reflection samples of the reflection pulse
% TDR_results.ERL               reported effective return loss
%
db = @(x) 20*log10(abs(x));
rms =@(x) norm(x)/sqrt(length(x));
if isfield(OP,'TDR_duration')
    TDR_duration=OP.TDR_duration; % approximate transit time multipler
else
    TDR_duration=5;
end
if ~isfield(OP,'DISPLAY_WINDOW')
    OP.DISPLAY_WINDOW=1; % approximate transit time multipler
end
f=sdd.Frequencies;
TDR_results.f=f;
% OP.Zt_adj=2;
if param.FLAG.S2P == 0
    
    % re-normalize reference of s-parameterss: this seems correct for a s4p input file
    TDR_RL =@(Zin,Zout,s11,s12,s21,s22)(Zin.^2.*s11+Zin.^2.*s22+Zout.^2.*s11+Zout.^2.*s22+Zin.^2-Zout.^2+Zin.*Zout.*s11.*2.0-Zin.*Zout.*s22.*2.0+Zin.^2.*s11.*s22-Zin.^2.*s12.*s21-Zout.^2.*s11.*s22+Zout.^2.*s12.*s21)./(Zin.*Zout.*2.0+Zin.^2.*s11+Zin.^2.*s22-Zout.^2.*s11-Zout.^2.*s22+Zin.^2+Zout.^2+Zin.^2.*s11.*s22-Zin.^2.*s12.*s21+Zout.^2.*s11.*s22-Zout.^2.*s12.*s21-Zin.*Zout.*s11.*s22.*2.0+Zin.*Zout.*s12.*s21.*2.0);
    
    if param.RL_sel==1, other_port=2;end
    if param.RL_sel==2, other_port=1;end
    for i = 1:length(sdd.Frequencies)
        if size(sdd.Parameters,2) ==1 % for s2p files
            RL(i)=TDR_RL(sdd.Impedance,2*ZT,sdd.Parameters(param.RL_sel,param.RL_sel,i) ,1, 1 ,sdd.Parameters(param.RL_sel,param.RL_sel,i) );
        else
            RL(i)=TDR_RL(sdd.Impedance,2*ZT,sdd.Parameters(param.RL_sel,param.RL_sel,i) ,sdd.Parameters( 1,2 ,i), sdd.Parameters( 2,1 ,i),sdd.Parameters( other_port,other_port ,i) );
        end
    end
    % elseif OP.Zt_adj ==2 % only adjust z_t drive impedance
    %     RL=squeeze((sdd.Parameters(param.RL_sel,param.RL_sel,:)));
    %     Z_t=ZT;
    %     zref=sdd.Impedance/2;
    %     if Z_t > zref
    %         radjust=  (zref-Z_t);
    %         S11adjust= radjust./(radjust + 2*zref);
    %         RL=RL  +S11adjust;
    %     elseif Z_t < zref
    %         rpad=-Z_t*zref/(Z_t-zref);
    %         S11adjust=zref/(rpad*(zref/rpad + 2));
    %         RL=RL   + S11adjust;
    %     else
    %         RL=RL;
    %     end
else
    for i = 1:length(sdd.Frequencies)
        rho= (2*ZT - sdd.Impedance) / (2*ZT + sdd.Impedance);
        interim =  sqrt(1-abs(rho)^2) * (1 - rho) / abs(1-rho);
        RL(i) = interim \ (sdd.Parameters(param.RL_sel,param.RL_sel,i) - rho) / ...
            (1 - rho *sdd.Parameters(param.RL_sel,param.RL_sel,i) ) * interim;
    end
end

% end
RL=squeeze(RL);
f9=f/1e9;
tr=param.TR_TDR;
TDR_results.delay=500e-12 ;
% determine max time from thue
% if sdd.NumPorts==1
%     try
%         maxtime=OP.N*param.ui;
%     catch
%         maxtime=2e-9;
%     end
%     pix=1;
% else
%     [ fir4del, tu] =get_RAW_FSIR(squeeze(sdd.Parameters(2,1,:)),f,OP,param);
%     pix=find(fir4del==max(fir4del),1);
%     maxtime=tu(pix)*TDR_duration+TDR_results.delay;
%     if maxtime > tu(end); maxtime=tu(end);end
% endS

try
    maxtime=OP.N*param.ui;
catch
    maxtime=2e-9;
end
if OP.N==0
    if sdd.NumPorts==1
        fprintf('<strong> Warning for s2p files N must not be zero<\strong> ');
    else
        [ fir4del, tu] =get_RAW_FIR(squeeze(sdd.Parameters(2,1,:)),f,OP,param);
        pix=find(fir4del==max(fir4del),1);
        maxtime=tu(pix)*TDR_duration+TDR_results.delay;
        if maxtime > tu(end); maxtime=tu(end);end
    end
end


% add delay 500 ps for TDR and 3 times Gaussnan transtion time
% (makes gausian edge somewhat causal)
H_t = exp( -2*(pi*f9*(tr)/1.6832).^2 ).*exp(-(1j)*2*pi*f9*TDR_results.delay/1e-9).*exp(-1j*2*pi*f9*tr*3);
if ~isfield(OP,'cb_Guassian')
    Use_gaussian=1;
else
    Use_gaussian=OP.cb_Guassian;
end
if Use_gaussian
    if iscolumn(H_t), H_t=H_t.'; end
    RLf=RL(:).'.*H_t;
else % add extra 3x tr delay for causality
    RLf=RL.*exp(-(1j)*2*pi*f9*TDR_results.delay/1e-9);
end

%Bessesl-Thomson turned off here (3rd input=0)
H_bt=Bessel_Thomson_Filter(param,f,0);

if isfield(OP,'TDR_Butterworth')
    H_bw=Butterworth_Filter(param,f,OP.TDR_Butterworth);
else
    H_bw=ones(1,length(f));
end


if param.Tukey_Window ~= 0
    H_tw= Tukey_Window(f,param);
else
    H_tw=ones(1,length(f));
end


if iscolumn(H_tw), H_tw=H_tw.';end
if iscolumn(H_bt), H_bt=H_bt.';end
if iscolumn(H_bw), H_bw=H_bw.';end
if iscolumn(RLf), RLf=RLf.';end

TDR_results.Rx_filter=H_bt.*H_bw.*H_tw;
RLf=RLf.*TDR_results.Rx_filter;
TDR_results.tx_filter=H_t;


[IR, t, causality_correction_dB, truncation_dB] = ...
    s21_to_impulse_DC(RLf,  sdd.Frequencies(:), param.sample_dt,OP);


%
% param.tfx =4.2e-10; % need to put in xls file and comment this out
tfx=param.tfx(np); % use fixture delay for port (np)

%  IR(abs(IR)<=OP.impulse_response_truncation_threshold)=0; % noise filter

t = t-TDR_results.delay;
tend=find(t>=maxtime+tfx,1); % n starts at tfx
if isempty(tend), tend=length(t); end
IR=IR(1:tend);
t=t(1:tend);
if isempty(tend), tend=length(t); end
tstart=find(t>=tr*1e-9,1); % account for Gaussian precursor
if isempty(tstart), tstart=1;end
if isempty(tend) || tstart >= tend
    if isempty(tend) || tstart >= tend
        %         warndlg('TDR compuation not valid, try decreasing truncation tolerance, increasing samples, or adding a transmisson line','WrnTDR');
    end
    tend=length(t);
    tstart=1;
end
OP.cb_step=0; % step is a basically a cos^2 form -pi/4 to 0. not activated.
ch=get_StepR(IR(tstart:tend),param,OP.cb_step,ZT);
TDR_results.tdr= ch.ZSR;
TDR_results.t = t(tstart:tend);

PTDR=get_PulseR(IR(tstart:tend),param,OP.cb_step,ZT);
if OP.TDR || OP.PTDR % determin average impededance with
    try
        tfstart=find(t>=3*tr*1e-9,1);
        x=squeeze(TDR_results.t(tfstart:end));TDR_results.x=TDR_results.tdr(:);
        y=squeeze(TDR_results.tdr(tfstart:end));TDR_results.y=TDR_results.t(:);
        w= exp(-(x-x(1))/OP.T_k ) ; % weighting function
        TDR_results.avgZport=mean(y.*w.')/mean(w.');
    catch
        TDR_results.avgZport=0;
        fit=zeros(1,1);
        p=[0 0 0 0  ];
    end
    TDR_results.RL=RL;
end
if OP.PTDR
    %     param.N_bx=param.ndfe;
    RL_equiv=-inf;
    L=param.levels;
    BinSize=OP.BinSize;
    %     param.specBER=1e-5;
    if OP.DISPLAY_WINDOW
        hwaitbar=waitbar(0);
    else
        fprintf('Worst ERL searching');
    end
    % adjust PTDR for NDFE
    % ---------------------- 2.7 code
    %     ntx=find(TDR_results.t >= tfx,1,'first');
    %     %     gatestartt=TDR_results.t(ntx);
    %     %     gatestartV=PTDR.pulse(ntx);
    %     ndfex=find(TDR_results.t > (param.N_bx+1)*param.ui+tfx,1,'first');
    %     tk=param.ui*1*(param.N_bx+1)+tfx;
    % -------------------
    % [ahealey 09/06/2019] Need to account for the 3*TR_TDR delay included in the rise
    % time filter.
    % ntx=find(TDR_results.t >= tfx,1,'first');
    ntx=find(TDR_results.t >= tfx+3*tr*1e-9,1,'first');
    %     gatestartt=TDR_results.t(ntx);
    %     gatestartV=PTDR.pulse(ntx);
    % ndfex=find(TDR_results.t > (param.N_bx+1)*param.ui+tfx,1,'first');
    % tk=param.ui*1*(param.N_bx+1)+tfx;
    ndfex=find(TDR_results.t > (param.N_bx+1)*param.ui+tfx+3*tr*1e-9,1,'first');
    tk=param.ui*1*(param.N_bx+1)+tfx+3*tr*1e-9;
    % [ahealey] End of modifications.
    if  isempty(ndfex), ndfex=length(TDR_results.t); end
    PTDR.pulse_orig=PTDR.pulse;
    
    switch param.Grr
        case 0 % pre .3cd release
            fctrx(1:length(PTDR.pulse_orig))=(1+param.rho_x)*param.rho_x;
        case 1 % .3cd release
            fctrx(1:length(PTDR.pulse_orig))=1;
        case 2 % .3ck working
            fctrx(1:length(PTDR.pulse_orig))=1;
    end
    Gloss(1:length(TDR_results.t))=1;
    Grr(1:length(TDR_results.t))=1;
    fctrx(1:ntx)=0; % moved out of loop for beta x and n_bx
    
    for ii=ntx:ndfex
        % adjust for near end loss
        if param.N_bx>0 && param.beta_x~=0;
            Gloss(ii)= 10.^(param.beta_x*(TDR_results.t(ii)-tk)/20);
        else
            Gloss(ii)=1;
        end
        % ---------------------- 2.7 code
        %         x=(TDR_results.t(ii)-tfx)/param.ui;
        % ----------------------
        % [ahealey9/06/2019] Need to account for the 3*TR_TDR delay included in the
        % rise time filter.
        % x=(TDR_results.t(ii)-tfx)/param.ui;
        x=(TDR_results.t(ii)-tfx-3*tr*1e-9)/param.ui;
        % determine how much of the return loss to use base on expected
        % missing reflections
        switch param.Grr
            case 0 % pre .3cd release
                Grr(ii)= (1+param.rho_x)*param.rho_x*exp(-(x-1*param.N_bx-1).^2/(1+param.N_bx)^2);
            case 1 % .3cd release
                Grr(ii)= (1+param.rho_x)*param.rho_x*exp(-(x-1*param.N_bx-1).^2/(1+param.N_bx)^2);
            case 2 % .3ck working
                Grr(ii)= param.rho_x ;
        end
        fctrx(ii)=Gloss(ii).*Grr(ii);
    end
    
    if isrow(fctrx), fctrx=fctrx(:);end
    PTDR.pulse=PTDR.pulse.*fctrx;
    if 0
        figure(10101+param.RL_sel);set(gcf,'Tag','COM');
        s1=subplot(2,1,1);
        plot((TDR_results.t(ntx:end)- TDR_results.t(ntx))/param.ui,fctrx(ntx:end),'disp','G_l_o_s_s*G_r_r');
        hold on
        plot((TDR_results.t(ntx:end)- TDR_results.t(ntx))/param.ui,Gloss(ntx:end),'disp','G_l_o_s_s');
        plot((TDR_results.t(ntx:end)- TDR_results.t(ntx))/param.ui,Grr(ntx:end),'disp','G_r_r');
        grid on
        ylim([ 0 1.2])
        s2=subplot(2,1,2);
        plot((TDR_results.t(ntx:end)- TDR_results.t(ntx))/param.ui,PTDR.pulse(ntx:end)./fctrx(ntx:end),'disp','PTDR');
        grid on
        linkaxes([s1,s2],'x')
        xlabel 'UI'
        xlim ([ 1 200])
    end
    
    FAST_NOISE_CONV=0;
    ERLRMS=rms(PTDR.pulse);
    for ki=1:param.samples_per_ui
        progress = ki/param.samples_per_ui;
        if OP.DISPLAY_WINDOW
            waitbar(progress, hwaitbar, 'ERL computing'); figure(hwaitbar); drawnow;
        else
            if ~mod(progress*100,1), fprintf('%i%% ', progress*100 );end
        end
        tps=PTDR.pulse(ki:param.samples_per_ui:end);
        if OP.RL_norm_test
            rl_fom=(norm(tps));
        else
            testpdf=get_pdf_from_sampled_signal( tps,L, BinSize*10, FAST_NOISE_CONV );
            cdf_test=cumsum(testpdf.y);
            rl_test=(-testpdf.x(find(cdf_test>=param.specBER,1,'first')));
            rl_fom=rl_test;
        end
        if rl_fom > RL_equiv
            RL_equiv=rl_fom;
            best_ki=ki;
        end
        if ~OP.RL_norm_test
            best_erl=rl_test;
            best_pdf=testpdf;
            best_cdf=cdf_test;
        end
        
    end
    if OP.RL_norm_test
        tps=PTDR.pulse(best_ki:param.samples_per_ui:end);
        testpdf=get_pdf_from_sampled_signal( tps,L, BinSize*10, FAST_NOISE_CONV );
        cdf_test=cumsum(testpdf.y);
        best_erl=(-testpdf.x(find(cdf_test>=param.specBER,1,'first')));
    end
    
    fprintf('\n');
    try
        close(hwaitbar)
    catch
    end
    if ~exist('best_ki','var'),best_ki=1;end
    TDR_results.ptdr_RL=PTDR.pulse;  % reflection waveform from the pulse
    TDR_results.WC_ptdr_samples_t=TDR_results.t(best_ki:param.samples_per_ui:end);
    TDR_results.WC_ptdr_samples=PTDR.pulse(best_ki:param.samples_per_ui:end);
    TDR_results.ERL=-db(best_erl);
    TDR_results.ERLRMS=-db(ERLRMS);
    
end


% end get TDR
%%


