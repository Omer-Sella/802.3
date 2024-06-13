function RILN_TD_struct= get_RILN_cmp_td(sdd21,RIL_struct,faxis_f2,OP,param,A_T)
% Complex reflection and re-reflection noise using the concept of zero'ing
% out of reflections
% sdd21 us a complex insertion loss
% RIL_struct is the output of capture_RIL_RILN()
% faxix_f2 needs to be at least to fb
% return reflections RILN_TD_struct.FOM based on time domain PR fit from pulse peak
% still need to settle on voltage scaling.
% maybe db(peak/Rss
db = @(x) 20*log10(abs(x));
fprintf('computing TD_RILN (dB) ...');

OP.interp_sparam_mag= 'trend_to_DC';
OP.interp_sparam_phase= 'interp_to_DC';
% OP.interp_sparam_mag= 'linear_trend_to_DC';
% OP.interp_sparam_phase= 'extrap_cubic_to_dc_linear_to_inf';

sdd21=squeeze(sdd21);
if iscolumn(sdd21)
    sdd21=sdd21.';
end
RIL=squeeze(RIL_struct.RIL);
if  iscolumn(RIL)
    RIL=RIL.';
end
rho_port1=squeeze(RIL_struct.rho_port1);
if  iscolumn(rho_port1)
    rho_port1=rho_port1.';
end
rho_port2=squeeze(RIL_struct.rho_port2);
if  iscolumn(rho_port2)
    rho_port2=rho_port2.';
end
RIL_f=squeeze(RIL_struct.freq);
if  iscolumn(RIL_f)
    RIL_f=RIL_f.';
end

%---start. Calculate the reflection and re-reflection noise
number_of_echos= 1e3;
fmin= 1e9;%<-------------
port2_reflection_rereflection_noise= zeros(1, length(RIL));
port1_reflection_rereflection_noise= zeros(1, length(RIL));
for m= 1:number_of_echos
    port2_reflection_rereflection_noise= port2_reflection_rereflection_noise+ abs(RIL).*(RIL.^(2*m)).*(rho_port1.^m).*(rho_port2.^m).*(1+rho_port1).*(1+rho_port2);
    port1_reflection_rereflection_noise= port1_reflection_rereflection_noise+ abs(RIL).*(RIL.^(2*m-1)).*(rho_port1.^(m-1)).*(rho_port2.^m).*(1+rho_port1).*(1+rho_port1);
end

%-----start. In the case of reflections, observed is bad TD conversion and hence removing data before 1GHz
fmin_idx= find(RIL_f>= fmin, 1, 'first');
port2_reflection_rereflection_noise= port2_reflection_rereflection_noise(fmin_idx:end);
port1_reflection_rereflection_noise= port1_reflection_rereflection_noise(fmin_idx:end);
f_reflection_rereflection_noise= RIL_f(fmin_idx:end);
%-----end. In the case of reflections, observed is bad TD conversion and hence removing data before 1GHz

% clear RIL RIL_f rho_port1 rho_port2
% clear fmin m
%---end. Calculate the reflection and re-reflection noise

fmbg=[ones(length(faxis_f2),1).*transpose(sdd21)  transpose(sqrt(faxis_f2)).*transpose(sdd21)  transpose(faxis_f2).*transpose(sdd21) transpose(faxis_f2.^2).*transpose(sdd21) ];
warning('off','MATLAB:nearlySingularMatrix');
unwraplog=log(abs(sdd21))+1i*unwrap(angle(sdd21));
LGw=transpose(sdd21.*unwraplog);
alpha = ((fmbg'*fmbg)^-1)*fmbg'*LGw;
efit_C=(alpha(1)+alpha(2).*sqrt(faxis_f2)+alpha(3).*faxis_f2 +faxis_f2.^2.*alpha(4)   );
FIT=transpose(exp(transpose(efit_C)));
efit=db(abs(FIT));
ILN = db(sdd21)-efit;


OP.impulse_response_truncation_threshold =1e-7;

print_for_codereview=0;
if exist('OP','var')
    H_bt=Bessel_Thomson_Filter(param,faxis_f2,1);
    H_bw=Butterworth_Filter(param,faxis_f2,1);
    H_t = exp(-(pi*faxis_f2/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
    H_tw=Tukey_Window(faxis_f2,param);
    H_tw=ones(1,length(faxis_f2) );
    [RILN_TD_struct.REF.FIR, ...
        RILN_TD_struct.REF.t, ...
        RILN_TD_struct.REF.causality_correction_dB, ...
        RILN_TD_struct.REF.truncation_dB] = s21_to_impulse_DC(sdd21.*H_bw.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    RILN_TD_struct.REF.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.REF.FIR);
    
        
    [RILN_TD_struct.FIT.FIR, ...
        RILN_TD_struct.FIT.t, ...
        RILN_TD_struct.FIT.causality_correction_dB, ...
        RILN_TD_struct.FIT.truncation_dB] = s21_to_impulse_DC(FIT.*H_bw.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    RILN_TD_struct.FIT.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.FIT.FIR);
    
    
    H_bt=Bessel_Thomson_Filter(param,RIL_f,1);
    H_bw=Butterworth_Filter(param,RIL_f,1);
    H_t = exp(-(pi*RIL_f/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
    H_tw=Tukey_Window(RIL_f,param);
    H_tw=ones(1,length(RIL_f) );
    [RILN_TD_struct.RIL.FIR, ...
        RILN_TD_struct.RIL.t, ...
        RILN_TD_struct.RIL.causality_correction_dB, ...
        RILN_TD_struct.RIL.truncation_dB] = s21_to_impulse_DC(RIL.*H_bw.*H_t.*H_tw ,RIL_f, param.sample_dt, OP) ;
    RILN_TD_struct.RIL.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.RIL.FIR);
    
    
    %---start. Calculate the channel delay
    try
    [delay_sec, delay_idx]= calculate_delay_CausalityEnforcement(faxis_f2, sdd21, param, OP);
    catch
    end
    port2_reflection_rereflection_noise= port2_reflection_rereflection_noise.*exp(-1j*2*pi*f_reflection_rereflection_noise*delay_sec);
    clear delay_sec delay_idx
    %---end. Calculate the channel delay


    
    H_bt=Bessel_Thomson_Filter(param,f_reflection_rereflection_noise,1);
    H_bw=Butterworth_Filter(param,f_reflection_rereflection_noise,1);
    H_t = exp(-(pi*f_reflection_rereflection_noise/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
    H_tw=Tukey_Window(f_reflection_rereflection_noise,param);
    H_tw=ones(1,length(f_reflection_rereflection_noise) );
    [RILN_TD_struct.REF_noise.FIR, ...
        RILN_TD_struct.REF_noise.t, ...
        RILN_TD_struct.REF_noise.causality_correction_dB, ...
        RILN_TD_struct.REF_noise.truncation_dB] = s21_to_impulse_DC(port2_reflection_rereflection_noise.*H_bw.*H_t.*H_tw ,f_reflection_rereflection_noise, param.sample_dt, OP) ;
    RILN_TD_struct.REF_noise.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.REF_noise.FIR);
    
    ipeak=find(RILN_TD_struct.REF.PR==max(RILN_TD_struct.REF.PR),1,'first');
    %     NrangeUI=1000;
    %     range_end=min(min(ipeak+param.samples_per_ui*NrangeUI,length(RILN_TD_struct.FIT.PR)-param.samples_per_ui ),length(RILN_TD_struct.REF.PR)-param.samples_per_ui);
    range_end= min(length(RILN_TD_struct.REF.PR), length(RILN_TD_struct.REF_noise.PR));
    range=ipeak:range_end;
    RILN_TD_struct.ILN=RILN_TD_struct.REF_noise.PR(range);
    RILN_TD_struct.t=RILN_TD_struct.REF_noise.t(range);
    RILN_TD_struct.FOM=-inf;
    RILN_TD_struct.FOM_PDF=-inf;
    rms_fom=-inf;
    for im=1:param.samples_per_ui
        RILN_TD_struct.FOM=max(RILN_TD_struct.FOM, norm( RILN_TD_struct.ILN(im:param.samples_per_ui:end)));
        [ pdf ] = get_pdf_from_sampled_signal(  RILN_TD_struct.ILN(im:param.samples_per_ui:end), param.levels, OP.BinSize ,0);
        rms=sqrt(pdf.y*pdf.x(:).^2)*sqrt(2);
        cdf=pdf; cdf.y=cumsum(pdf.y);
        %         cursors = d_cpdf(OP.BinSize,param.a_thru*[-1 1], [1 1]/2);
        %         signal_and_isi_pdf = conv_fct(cursors, pdf);
        %         cdf=signal_and_isi_pdf; cdf.y=cumsum(signal_and_isi_pdf.y);
        if print_for_codereview  % remove once all checked out
            h=figure(190);set(gcf,'Tag','COM');
            semilogy(-cdf.x,cdf.y);
            %             xlim ([0,-cdf.x(1)])
            ylim([param.specBER 1]);title ('CDF of ILN')
            hold on
        end
        if rms>rms_fom
            rms_fom=rms;
            RILN_TD_struct.FOM_PDF= -cdf.x(find(cdf.y >= param.specBER, 1, 'first'));
            RILN_TD_struct.PDF=pdf;
        end
    end
    pdf_from_norm=normal_dist(RILN_TD_struct.FOM, 7 , OP.BinSize);
    RILN_TD_struct.SNR_ISI_FOM=db(RILN_TD_struct.FIT.PR(ipeak)/RILN_TD_struct.FOM);
    RILN_TD_struct.SNR_ISI_FOM_PDF=db(RILN_TD_struct.FIT.PR(ipeak)/RILN_TD_struct.FOM_PDF);
    %     fprintf('%g dB\n',RILN_TD_struct.SNR_ISI_FOM)
    fprintf('%g dB\n',RILN_TD_struct.SNR_ISI_FOM_PDF)
    if print_for_codereview % remove once all checked out
        figure(9000);set(gcf,'Tag','COM');
        plot(RILN_TD_struct.REF.t,RILN_TD_struct.REF.PR,'disp','ref')
        hold on
        plot(RILN_TD_struct.FIT.t,RILN_TD_struct.FIT.PR,'disp','fit')
        plot(RILN_TD_struct.RIL.t,RILN_TD_struct.RIL.PR,'disp','RILN')
        yyaxis right
        plot(RILN_TD_struct.t,RILN_TD_struct.ILN,'disp','td RILN')
        hold off
        fprintf('SNR ISI FOM rms = %g dB;   SNR ISI FOM PDF = %g dB\n',RILN_TD_struct.SNR_ISI_FOM,RILN_TD_struct.SNR_ISI_FOM_PDF)
        figure(9002);set(gcf,'Tag','COM');
        semilogy(RILN_TD_struct.PDF.x,RILN_TD_struct.PDF.y,'disp','actual PDF')
        hold on
        semilogy(pdf_from_norm.x,pdf_from_norm.y,'disp','PDF using Gaussian assumed PDF');
        ylim([param.specBER max([RILN_TD_struct.PDF.y pdf_from_norm.y])]);title ('Compare actual PDF to Gaussian')
        grid on
        legend('show')
    end
end
function result=get_StepR(ir,param,cb_step,ZT)
%ir = impulse response
%t_base=time array with equal time steps
%samp_UI = number of samples per UI for ir
%   result.SBR
% t for debug
t=(1/param.fb)/param.samples_per_ui*(0:length(ir)-1);

if cb_step
    Ag=1;
    dt=1/param.fb/param.samples_per_ui;
    edge_time=param.TR_TDR*1e-9;
    fedge=1/edge_time;
    tedge=0:dt:edge_time*2;
    %
    edge=Ag*(2*cos(2*pi*(tedge)*fedge/16-pi/4).^2-1);
    drive_pulse=[edge ones(1,param.samples_per_ui)];
    %pulse=filter(UI_ones,1,ir);
    
    pulse=filter(drive_pulse,1,ir);
else
    pulse=cumsum(ir);
end
TDR_response=(1+pulse)./(1-pulse)*ZT*2;
result.ZSR=TDR_response;
result.pulse=pulse;
