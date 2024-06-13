function chdata=Apply_EQ(param,fom_result,chdata,OP)

FB=param.fb;
FZ=param.CTLE_fz(fom_result.ctle);
FP1=param.CTLE_fp1(fom_result.ctle);
FP2=param.CTLE_fp2(fom_result.ctle);
GDC=param.ctle_gdc_values(fom_result.ctle);
if ~isempty(param.f_HP)
    FHP=param.f_HP(fom_result.best_G_high_pass);
end
if ~isempty(param.g_DC_HP_values)
    GDCHP=param.g_DC_HP_values(fom_result.best_G_high_pass);
end
if ~isempty(param.f_HP_Z)
    FHPZ=param.f_HP_Z(fom_result.ctle);
end
if ~isempty(param.f_HP_P)
    FHPP=param.f_HP_P(fom_result.ctle);
end
%Handle the scenario where the pulse response is not long enough to
%contain all DFE taps.  the SBR recorded in fom_result has the proper
%length
SBR_Len=length(fom_result.sbr);
if length(chdata(1).uneq_imp_response)<SBR_Len
    samples_added=SBR_Len-length(chdata(1).uneq_imp_response);
    chdata(1).uneq_imp_response(end+1:SBR_Len)=0;
    chdata(1).uneq_pulse_response(end+1:SBR_Len)=0;
    chdata(1).t(end+1:SBR_Len)=(1:samples_added)/param.fb/param.samples_per_ui+chdata(1).t(end);
end
for i=1:param.number_of_s4p_files
    % get quick PDF  results but only for THRU when in Rx calibration
    uneq_ir=chdata(i).uneq_imp_response;% includes packages, Hx, and Hr
    if OP.INCLUDE_CTLE==1
        switch param.CTLE_type
            case 'CL93'
                eq_ir = TD_CTLE(uneq_ir, FB, FZ, FP1, FP2, GDC, param.samples_per_ui);
            case 'CL120d'
                eq_ir = TD_CTLE(uneq_ir, FB, FZ, FP1, FP2, GDC, param.samples_per_ui);
                eq_ir = TD_CTLE(eq_ir, FB, FHP, FHP, 100e100 , GDCHP, param.samples_per_ui);
            case 'CL120e' % z has been adjusted  for gain
                eq_ir = TD_CTLE(uneq_ir, FB, FZ, FP1, FP2, GDC, param.samples_per_ui);
                eq_ir = TD_CTLE(eq_ir,FB, FHPZ, FHPP,1e99, 0, param.samples_per_ui);
        end
    else
        eq_ir=uneq_ir;
    end
    chdata(i).eq_imp_response=eq_ir;
    eq_pulse=filter(ones(1, param.samples_per_ui), 1, chdata(i).eq_imp_response);

    if isequal(chdata(i).type, 'FEXT') || isequal(chdata(i).type, 'THRU')
        eq_pulse = FFE( fom_result.txffe ,fom_result.cur-1 , param.samples_per_ui,  eq_pulse );
    end
    % chdata(i).ctle_imp_response
    if OP.RxFFE
        if isequal(upper(OP.FFE_OPT_METHOD),'MMSE')
            chdata(i).ctle_imp_response = FFE( fom_result.RxFFE ,fom_result.cur-1 , param.samples_per_ui,  eq_ir );
        end
        [ eq_pulse, C]=force(eq_pulse,param,OP,fom_result.t_s,fom_result.RxFFE);
    end
    chdata(i).eq_pulse_response=eq_pulse;% includes packages, Hf, and Hr, Ht, and Hffe(from tx)
end