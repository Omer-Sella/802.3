function output_args=Output_Arg_Fill(output_args,sigma_bn,Noise_Struct,COM_SNR_Struct,param,chdata,fom_result,OP)

%not all output_args are filled here but most are

switch lower(OP.TDECQ)
    case { false 'none' } % should be the default
        output_args.VMA=[];
    case 'vma'
        est_vma=vma(fom_result.sbr,param.samples_per_ui);
        output_args.VMA=est_vma.VMA;
    otherwise
        error('%s not recognized for Histogram_Window_Weigh this feature is limited',OP.TDECQ)
end

fileset_str=str2csv({chdata.base});
output_args.file_names=sprintf('"%s"', fileset_str);
% [ahealey] Echo the termination parameters in the output arguments..
for odt_param = {'R_diepad', 'C_diepad', 'L_comp', 'C_bump'}
    output_args.(odt_param{:}) = param.(odt_param{:});
end
% [ahealey] End of modifications.
for pkg_params = {'levels', 'Pkg_len_TX', 'Pkg_len_NEXT', 'Pkg_len_FEXT', 'Pkg_len_RX','R_diepad','pkg_Z_c','C_v'}
    output_args.(pkg_params{:})= param.(pkg_params{:});
end
output_args.baud_rate_GHz=param.fb/1e9;
output_args.f_Nyquist_GHz = param.fb/2e9;
output_args.BER=param.specBER;
output_args.FOM = fom_result.FOM;
output_args.sigma_N=Noise_Struct.sigma_N;
output_args.DFE4_RSS=norm(fom_result.DFE_taps(4:end));
output_args.DFE2_RSS=norm(fom_result.DFE_taps(2:end));
output_args.tail_RSS=fom_result.tail_RSS;
output_args.channel_operating_margin_dB=COM_SNR_Struct.COM;
output_args.available_signal_after_eq_mV=1000*COM_SNR_Struct.A_s;
output_args.peak_uneq_pulse_mV=1000*max(abs(chdata(1).uneq_pulse_response));
try
    output_args.uneq_FIR_peak_time=chdata(1).t(chdata(1).uneq_imp_response==max(chdata(1).uneq_imp_response));
catch
    output_args.uneq_FIR_peak_time=[];
end
output_args.steady_state_voltage_mV = 1000*fom_result.A_f; % RIM 7/03/2019 use peak from optimize_FOM
its=find(chdata(1).eq_pulse_response>=max(chdata(1).eq_pulse_response),1,'first');
isumend=min(its+param.N_v*param.samples_per_ui,length(chdata(1).eq_pulse_response));
output_args.steady_state_voltage_weq_mV = 1000*sum(chdata(1).eq_pulse_response(1:isumend) )/param.samples_per_ui;

if OP.RX_CALIBRATION== 1
    output_args.sigma_bn=sigma_bn;
else
    output_args.sigma_bn=[];
end
output_args.Peak_ISI_XTK_and_Noise_interference_at_BER_mV=1000*COM_SNR_Struct.A_ni;
output_args.peak_ISI_XTK_interference_at_BER_mV=1000*Noise_Struct.peak_interference_at_BER;
output_args.peak_ISI_interference_at_BER_mV=1000*Noise_Struct.thru_peak_interference_at_BER;
output_args.equivalent_ICI_sigma_assuming_PDF_is_Gaussian_mV=Noise_Struct.sci_sigma*1000;

if OP.RX_CALIBRATION == 0
    output_args.peak_MDXTK_interference_at_BER_mV=1000*Noise_Struct.crosstalk_peak_interference_at_BER;
    output_args.peak_MDNEXT_interference_at_BER_mV=1000*Noise_Struct.MDNEXT_peak_interference;
    output_args.peak_MDFEXT_interference_at_BER_mV=1000*Noise_Struct.MDFEXT_peak_interference;
else
    output_args.peak_MDXTK_interference_at_BER_mV=[];
    output_args.peak_MDNEXT_interference_at_BER_mV=[];
    output_args.peak_MDFEXT_interference_at_BER_mV=[];
end
%output_args.ICN_mV=ICN*1000;
%             output_args.ICN_test_mV=ICN_test*1000;
xtk=param.num_next+param.num_fext;
if xtk>0 && OP.RX_CALIBRATION ==0 && OP.TDMODE==0
    %output_args.MDNEXT_ICN_92_46_mV=MDNEXT_ICN*1000;
    %output_args.MDFEXT_ICN_92_47_mV=MDFEXT_ICN*1000;
    output_args.equivalent_ICN_assuming_Gaussian_PDF_mV=Noise_Struct.cci_sigma*1000;
else
    output_args.MDNEXT_ICN_92_46_mV=0;
    output_args.MDFEXT_ICN_92_47_mV=0;
    output_args.equivalent_ICN_assuming_PDF_is_Gaussian_mV=0;
end
%output_args.SNR_ISI_XTK_normalized_1_sigma=20*log10(A_s/(peak_interference_at_BER/qfuncinv(param.specBER))); modified by Yasuo Hidaka, 8/7/17
if 1
    output_args.SNR_ISI_XTK_normalized_1_sigma=20*log10(COM_SNR_Struct.A_s/(Noise_Struct.peak_interference_at_BER/sqrt(2)/erfcinv(2*param.specBER)));
    output_args.SNR_ISI_est=fom_result.SNR_ISI;
    output_args.Pmax_by_Vf_est=fom_result.Pmax_by_Vf;
    output_args.Tr_measured_from_step_ps=fom_result.Tr_measured_from_step/1e-12;
end


switch param.CTLE_type
    case 'CL93'
        output_args.CTLE_zero_poles=[param.CTLE_fz(fom_result.ctle) param.CTLE_fp2(fom_result.ctle) param.CTLE_fp1(fom_result.ctle)];
        output_args.CTLE_DC_gain_dB=param.ctle_gdc_values(fom_result.ctle);
        output_args.g_DC_HP=[];
        output_args.HP_poles_zero=[];
    case 'CL120d'
        output_args.CTLE_zero_poles=[param.CTLE_fz(fom_result.ctle) param.CTLE_fp2(fom_result.ctle) param.CTLE_fp1(fom_result.ctle)];
        output_args.CTLE_DC_gain_dB=param.ctle_gdc_values(fom_result.ctle);
        output_args.g_DC_HP=param.g_DC_HP_values(fom_result.best_G_high_pass);
        output_args.HP_poles_zero=param.f_HP(fom_result.best_G_high_pass);
    case 'CL120e'
        output_args.CTLE_zero_poles=[param.CTLE_fz(fom_result.ctle) param.f_HP_Z(fom_result.ctle) param.CTLE_fp2(fom_result.ctle) param.CTLE_fp1(fom_result.ctle) param.f_HP_P(fom_result.ctle)];
        output_args.CTLE_DC_gain_dB=param.ctle_gdc_values(fom_result.ctle);
        output_args.g_DC_HP=[];
        output_args.HP_poles_zero=[];
end
output_args.TXLE_taps=fom_result.txffe;
if  length(output_args.TXLE_taps) >= 3
    output_args.Pre2Pmax = -output_args.TXLE_taps(end-2)/output_args.TXLE_taps(end-1);
else
    output_args.Pre2Pmax=[];
end
output_args.DFE_taps=fom_result.DFE_taps;
if param.Floating_DFE ||  param.Floating_RXFFE
    output_args.floating_tap_locations=fom_result.floating_tap_locations;
else
    output_args.floating_tap_locations=[];
end

if OP.RxFFE
    output_args.RxFFE=fom_result.RxFFE;
    output_args.RxFFEgain=param.current_ffegain;
else % Yasou Hidaka 11/20/2018 help to align csv file columns
    output_args.RxFFE=[];
    output_args.RxFFEgain=[];
end

output_args.itick=fom_result.itick;

% Calculation of error propagation and burst probability
if OP.nburst>0
    [p_burst,p_error_propagation]=Burst_Probability_Calc(COM_SNR_Struct,fom_result.DFE_taps,param,OP);
    output_args.error_propagation_probability = p_error_propagation;
    output_args.burst_probabilities = p_burst;
else
    output_args.error_propagation_probability = [];
    output_args.burst_probabilities = [];
end


%begin yasuo patch 12/11/2018
% collect sigma values to report
% pdf2sgm() is a function to calculate sigma value from PDF
% It is added at the end of this file code.
% I am not sure if an equivalent function already exists.
output_args.sgm_Ani__isi_xt_noise = pdf2sgm(COM_SNR_Struct.combined_interference_and_noise_pdf);
output_args.sgm_isi_xt = pdf2sgm(Noise_Struct.isi_and_xtalk_pdf);
output_args.sgm_noise__gaussian_noise_p_DD  = pdf2sgm(Noise_Struct.noise_pdf);
output_args.sgm_p_DD = pdf2sgm(Noise_Struct.p_DD);
output_args.sgm_gaussian_noise = pdf2sgm(Noise_Struct.gaussian_noise_pdf);
output_args.sgm_G = Noise_Struct.sigma_G;
output_args.sgm_rjit = Noise_Struct.sigma_rjit;
output_args.sgm_N = Noise_Struct.sigma_N;
output_args.sgm_TX = Noise_Struct.sigma_TX;
output_args.sgm_isi = pdf2sgm(Noise_Struct.sci_pdf);
if OP.RX_CALIBRATION == 0
    output_args.sgm_xt  = pdf2sgm(Noise_Struct.cci_pdf);
else
    output_args.sgm_xt=[];
end
% end yasuo patch

% r259 putting COM, VEO and loss last in report
%         output_args.VEO_normalized = (A_s-A_ni)/A_s;
output_args.VEC_dB = COM_SNR_Struct.VEC_dB;
output_args.VEO_mV = COM_SNR_Struct.VEO_mV;
if OP.RX_CALIBRATION ==0 && OP.EW == 1
    output_args.EW_UI_est=COM_SNR_Struct.EW_UI;
    output_args.eye_contour=COM_SNR_Struct.eye_contour;
    output_args.VEO_window_mUI= param.T_O;
else
    output_args.EW_UI_est=[];
    output_args.eye_contour=[];
    output_args.VEO_window_mUI= [];
end

if sum(param.AC_CM_RMS) ~= 0
    output_args.sigma_ACCM_at_tp0_mV=chdata(1).sigma_ACCM_at_tp0*1000;
    fprintf(' AC RMS at TP0 = %.3g mV \n',output_args.sigma_ACCM_at_tp0_mV)
    output_args.sigma_AC_CCM_at_rxpkg_output_mV=chdata(1).CD_CM_RMS*1000; %
else
    output_args.sigma_ACCM_at_tp0_mV=[];
    output_args.sigma_AC_CCM_at_rxpkg_output_mV=[];
end
if OP.MLSE
    output_args.COM_orig=COM_SNR_Struct.COM_orig;
    output_args.VEC_dB_orig=COM_SNR_Struct.VEC_dB_orig;
end
%
output_args.COM_dB=COM_SNR_Struct.COM;
% end yasuo patch
% begin yasuo patch 3/18/2019
output_args.DER_thresh = COM_SNR_Struct.threshold_DER;
% end yasuo patch