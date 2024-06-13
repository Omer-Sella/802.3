function [PDF,CDF,NS]=Create_Noise_PDF(A_s,param,fom_result,chdata,OP,sigma_bn,PSD_results)

%This block was originally in main COM function but was moved here for
%cleanup.  It returns the combined interference and noise PDF & CDF as well
%as a structure "NS" that contains all the noise parameters that are used
%in other places in COM

if OP.RX_CALIBRATION
    ctle_gain2 = (10.^(param.ctle_gdc_values(fom_result.ctle)/20) + 1i*chdata(2).faxis/param.CTLE_fz(fom_result.ctle)) ./ ...
        ((1+1i*chdata(2).faxis/param.CTLE_fp1(fom_result.ctle)).*(1+1i*chdata(2).faxis/param.CTLE_fp2(fom_result.ctle)));
    switch param.CTLE_type
        case 'CL93'
            H_low2=1;
        case 'CL120d' % this clause uses two gain indexes
            H_low2=(10.^(param.g_DC_HP_values(fom_result.best_G_high_pass)/20) +  1i*chdata(2).faxis/param.f_HP(fom_result.best_G_high_pass))./(1 + 1i*chdata(2).faxis/param.f_HP(fom_result.best_G_high_pass));
        case 'CL120e' % Z1 has been adjusted
            H_low2=(1 +  1i*chdata(2).faxis/f_HP_P(fom_result.ctle))./(1 + 1i*chdata(2).faxis/f_HP_Z(fom_result.ctle));
    end
    H_ctf2=H_low2.*ctle_gain2;
    [ sigma_ne, NS.sigma_hp] = get_sigma_noise( H_ctf2,  param, chdata, sigma_bn );
else
    sigma_ne=0;
end

NS.sigma_N = fom_result.sigma_N; % eta zero noise
if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
    if ~OP.SNR_TXwC0
        % Equation 93A-30 %% h0(ts0)= A_s*R_lm/(L-1)
        NS.sigma_TX = (param.levels-1)*A_s/param.R_LM*10^(-param.SNR_TX/20); % SNR_Tx and RLM
    else
        NS.sigma_TX = (param.levels-1)*A_s/fom_result.txffe(fom_result.cur) /param.R_LM*10^(-param.SNR_TX/20); % SNR_Tx from Adee
    end
else
    NS.sigma_TX =PSD_results.tn_rms;
end
% Equation 93A-41 %%
NS.sigma_G  = norm([param.sigma_RJ*param.sigma_X*norm(fom_result.h_J), NS.sigma_N, NS.sigma_TX]);
NS.sigma_rjit= param.sigma_RJ*param.sigma_X*norm(fom_result.h_J);

% Equation 93A-42 %%
% number of sigmas needed depends on the required BER.
if param.Noise_Crest_Factor == 0
    NS.ber_q = sqrt(2)*erfcinv(2*param.specBER);
else
    NS.ber_q=param.Noise_Crest_Factor;
end
NS.gaussian_noise_pdf = normal_dist(NS.sigma_G, NS.ber_q, param.delta_y);
% enable overriding the Q factor of the BBN instrument.
if OP.force_BBN_Q_factor
    NS.ne_noise_pdf = normal_dist(sigma_ne, OP.BBN_Q_factor, param.delta_y);
else
    NS.ne_noise_pdf = normal_dist(sigma_ne, NS.ber_q, param.delta_y);
end
NS.gaussian_noise_pdf = conv_fct(NS.gaussian_noise_pdf, NS.ne_noise_pdf);

% p_DD is computed using the procedure defined in 93A.1.7.1 with h(n)=A_DD*h_J(n)
NS.p_DD = get_pdf_from_sampled_signal(param.A_DD*fom_result.h_J, param.levels, param.delta_y);

% Equation 93A-43
NS.noise_pdf=conv_fct(NS.gaussian_noise_pdf, NS.p_DD);

gaussian_rjitt_pdf = normal_dist(NS.sigma_rjit, NS.ber_q, param.delta_y);
NS.jitt_pdf=conv_fct(gaussian_rjitt_pdf, NS.p_DD);

% Implementation of 93A.1.7.3 combination procedure
%  (effectively Equation 93A-44) %%

% Self-Channel Interference is thru residual result
NS.sci_pdf = chdata(1).pdfr;
sci_mxi=find(cumsum(NS.sci_pdf.y)>=param.specBER, 1, 'first');
NS.thru_peak_interference_at_BER=abs(NS.sci_pdf.x(sci_mxi));
sci_msi=find(cumsum(NS.sci_pdf.y)>=param.specBER, 1, 'first');
NS.sci_sigma=abs(NS.sci_pdf.x(sci_msi)/(erfcinv(2*param.specBER)*sqrt(2)));
if OP.RX_CALIBRATION ==0
    % Co-Channel Interference PDFs (for information only):
    % initialize to deltas
    MDNEXT_cci_pdf = d_cpdf(param.delta_y, 0, 1);
    MDFEXT_cci_pdf = d_cpdf(param.delta_y, 0, 1);
    % serially convolve FEXT/NEXT PDFs
    for k=2:param.number_of_s4p_files
        if isequal(chdata(k).type, 'NEXT')
            MDNEXT_cci_pdf = conv_fct(MDNEXT_cci_pdf, chdata(k).pdfr);
        else % ... must be FEXT
            MDFEXT_cci_pdf = conv_fct(MDFEXT_cci_pdf, chdata(k).pdfr);
        end
    end
    
    % find "peaks" of MDNEXT/MDFEXT for reporting
    mdnxi=find(cumsum(MDNEXT_cci_pdf.y)>=param.specBER, 1, 'first');
    NS.MDNEXT_peak_interference=abs(MDNEXT_cci_pdf.x(mdnxi));
    mdfxi=find(cumsum(MDFEXT_cci_pdf.y)>=param.specBER, 1, 'first');
    NS.MDFEXT_peak_interference=abs(MDFEXT_cci_pdf.x(mdfxi));
    
    % Combined crosstalk effect
    NS.cci_pdf = conv_fct(MDFEXT_cci_pdf, MDNEXT_cci_pdf);
    cci_mxi=find(cumsum(NS.cci_pdf.y)>=param.specBER, 1, 'first');
    cci_msi=find(cumsum(NS.cci_pdf.y)>=param.specBER, 1, 'first');
    NS.cci_sigma=abs(NS.cci_pdf.x(cci_msi)/(erfcinv(2*param.specBER)*sqrt(2)));
    NS.crosstalk_peak_interference_at_BER=abs(NS.cci_pdf.x(cci_mxi));
    % combine cci and sci
    NS.isi_and_xtalk_pdf = conv_fct(NS.sci_pdf, NS.cci_pdf);
else
    % for calibration there is no cci
    NS.isi_and_xtalk_pdf=NS.sci_pdf;
end

mxi=find(cumsum(NS.isi_and_xtalk_pdf.y)>=param.specBER, 1, 'first');
NS.peak_interference_at_BER=abs(NS.isi_and_xtalk_pdf.x(mxi));


% Equation 93A-45
combined_interference_and_noise_pdf = conv_fct(NS.isi_and_xtalk_pdf, NS.noise_pdf);
PDF=combined_interference_and_noise_pdf;

% Equation 93A-37
combined_interference_and_noise_cdf=cumsum(combined_interference_and_noise_pdf.y);
CDF=combined_interference_and_noise_cdf;