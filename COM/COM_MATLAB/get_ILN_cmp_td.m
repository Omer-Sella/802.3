function [ILN, efit, TD_ILN]= get_ILN_cmp_td(sdd21,faxis_f2,OP,param,A_T)
% Complex IL fitting
% sdd21 us a complex insertion loss
% efit and ILN are in db
% faxix_f2 needs to be at least to fb
% return reflections TD_ILN.FOM based on time domain PR fit from pulse peak
% still need to settle on voltage scaling.
% maybe db(peak/Rss

OP.interp_sparam_mag= 'trend_to_DC';
OP.interp_sparam_phase= 'interp_to_DC';
% OP.interp_sparam_mag= 'linear_trend_to_DC';
% OP.interp_sparam_phase= 'extrap_cubic_to_dc_linear_to_inf';

print_for_codereview=0;
if ~exist('A_T','var')
    A_T=1;
end

db = @(x) 20*log10(abs(x));
sdd21=squeeze(sdd21);
if  iscolumn(sdd21)
    sdd21=sdd21.';
end
fmbg=[ones(length(faxis_f2),1).*transpose(sdd21)  transpose(sqrt(faxis_f2)).*transpose(sdd21)  transpose(faxis_f2).*transpose(sdd21) transpose(faxis_f2.^2).*transpose(sdd21) ];
warning('off','MATLAB:nearlySingularMatrix');
unwraplog=log(abs(sdd21))+1i*unwrap(angle(sdd21));
LGw=transpose(sdd21.*unwraplog);
alpha = ((fmbg'*fmbg)^-1)*fmbg'*LGw;
efit_C=(alpha(1)+alpha(2).*sqrt(faxis_f2)+alpha(3).*faxis_f2 +faxis_f2.^2.*alpha(4)   );
FIT=transpose(exp(transpose(efit_C)));
efit=db(abs(FIT));
ILN = db(sdd21)-efit;
% time domain
fprintf('computing TD_ILN (dB) ...')
if exist('OP','var')
%     OP.fraction_of_F_range_start_extrap_from=.95;
    OP.impulse_response_truncation_threshold =1e-7;
    
    H_bt=Bessel_Thomson_Filter(param,faxis_f2,1);
    H_bw=Butterworth_Filter(param,faxis_f2,1);
    H_t = exp(-(pi*faxis_f2/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
    H_tw=Tukey_Window(faxis_f2,param);
    H_tw=ones(1,length(faxis_f2) );
    
   [TD_ILN.REF.FIR, ...
        TD_ILN.REF.t, ...
        TD_ILN.REF.causality_correction_dB, ...
        TD_ILN.REF.truncation_dB] = s21_to_impulse_DC(sdd21.*H_bt.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    TD_ILN.REF.PR=filter(ones(1, param.samples_per_ui), 1, TD_ILN.REF.FIR);

    [TD_ILN.FIT.FIR, ...
        TD_ILN.FIT.t, ...
        TD_ILN.FIT.causality_correction_dB, ...
        TD_ILN.FIT.truncation_dB] = s21_to_impulse_DC(FIT.*H_bt.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    TD_ILN.FIT.PR=filter(ones(1, param.samples_per_ui), 1, TD_ILN.FIT.FIR);
    ipeak=find(TD_ILN.REF.PR==max(TD_ILN.REF.PR),1,'first');
    %     NrangeUI=1000;
    %     range_end=min(min(ipeak+param.samples_per_ui*NrangeUI,length(TD_ILN.FIT.PR)-param.samples_per_ui ),length(TD_ILN.REF.PR)-param.samples_per_ui);
    range_end= min(length(TD_ILN.REF.PR), length(TD_ILN.FIT.PR));
    range=ipeak:range_end;
    TD_ILN.ILN=TD_ILN.FIT.PR(range)-TD_ILN.REF.PR(range);
    TD_ILN.t=TD_ILN.FIT.t(range);
    TD_ILN.FOM=-inf;
    TD_ILN.FOM_PDF=-inf;
    rms_fom=-inf;
    for im=1:param.samples_per_ui
        TD_ILN.FOM=max(TD_ILN.FOM, norm( TD_ILN.ILN(im:param.samples_per_ui:end)));
        [ pdf ] = get_pdf_from_sampled_signal(  TD_ILN.ILN(im:param.samples_per_ui:end), param.levels, OP.BinSize ,0);
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
            TD_ILN.FOM_PDF= -cdf.x(find(cdf.y >= param.specBER, 1, 'first'));
            TD_ILN.PDF=pdf;
        end
    end
    pdf_from_norm=normal_dist(TD_ILN.FOM, 7 , OP.BinSize);
    TD_ILN.SNR_ISI_FOM=db(TD_ILN.FIT.PR(ipeak)/TD_ILN.FOM);
    TD_ILN.SNR_ISI_FOM_PDF=db(TD_ILN.FIT.PR(ipeak)/TD_ILN.FOM_PDF);
    %     fprintf('%g dB\n',TD_ILN.SNR_ISI_FOM)
    fprintf('%g dB\n',TD_ILN.SNR_ISI_FOM_PDF)
    if print_for_codereview % remove once all checked out
        figure(9000);set(gcf,'Tag','COM');
        plot(TD_ILN.t,TD_ILN.ILN,'disp','td iln')
        hold on
        plot(TD_ILN.FIT.t,TD_ILN.FIT.PR,'disp','fit')
        plot(TD_ILN.REF.t,TD_ILN.REF.PR,'disp','ref')
        hold off
        fprintf('SNR ISI FOM rms = %g dB;   SNR ISI FOM PDF = %g dB\n',TD_ILN.SNR_ISI_FOM,TD_ILN.SNR_ISI_FOM_PDF)
        figure(9002);set(gcf,'Tag','COM');
        semilogy(TD_ILN.PDF.x,TD_ILN.PDF.y,'disp','actual PDF')
        hold on
        semilogy(pdf_from_norm.x,pdf_from_norm.y,'disp','PDF using Gaussian assumed PDF');
        ylim([param.specBER max([TD_ILN.PDF.y pdf_from_norm.y])]);title ('Compare actual PDF to Gaussian')
        grid on
        legend('show')
    end
end
% display('got to end of get_ILN_cmp_td')