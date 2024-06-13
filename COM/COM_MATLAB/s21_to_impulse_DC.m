function [voltage, t_base, causality_correction_dB, truncation_dB] = ...
    s21_to_impulse_DC(IL, freq_array, time_step, OP)
% Creates a time-domain impulse response from frequency-domain IL data.
% IL does not need to have DC but a corresponding frequency array
% (freq_array) is required.
%
% Causality is imposed using the Alternating Projections Method. See also:
% Quatieri and Oppenheim, "Iterative Techniques for Minimum Phase Signal
% Reconstruction from Phase or Magnitude", IEEE Trans. ASSP-29, December
% 1981 (http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163714)

ILin=IL;
fmax=1/time_step/2;
freq_step=(freq_array(3)-freq_array(2))/1;
fout=0:1/round(fmax/freq_step)*fmax:fmax;
if all(IL==0)
    %response with all zeros is problematic.  set to all eps and avoid interp function
    IL=ones(1,length(fout))*eps;
else
    IL=interp_Sparam(ILin,freq_array,fout, OP.interp_sparam_mag, OP.interp_sparam_phase,OP);
    IL_nan = find(isnan(IL));
    for in=IL_nan
        IL(in)=IL(in-1);
    end
end
IL = IL(:);
% add padding for time steps
% IL_symmetric = [IL(1:end-1);0; flipud(conj(IL(2:end-1)))];
IL_symmetric = [real(IL(1)); IL(2:end-1); real(IL(end)); flipud(conj(IL(2:end-1)))];
impulse_response = real(ifft(IL_symmetric));
L = length(impulse_response);
t_base = (0:L-1)/(freq_step*L);

original_impulse_response=impulse_response;
% Correct non-causal effects frequently caused by extrapolation of IL
% Assumption: peak of impulse_response is in the first half, i.e. not anti-causal
abs_ir=abs(impulse_response);
a = find(abs_ir(1:L/2) > max(abs_ir(1:L/2))*OP.EC_PULSE_TOL);
start_ind = a(1);

err=inf;
while ~all(impulse_response==0)
    impulse_response(1:start_ind)=0;
    impulse_response(floor(L/2):end)=0;
    IL_modified=abs(IL_symmetric).*exp(1j*angle(fft(impulse_response)));
    ir_modified = real(ifft(IL_modified));
    delta = abs(impulse_response-ir_modified);
    
    err_prev = err;
    err=max(delta)/max(impulse_response);
    if err<OP.EC_REL_TOL || abs(err_prev-err)<OP.EC_DIFF_TOL
        break;
    end
    
    impulse_response=ir_modified;
end

causality_correction_dB=20*log10(norm(impulse_response-original_impulse_response)/norm(impulse_response));

if ~OP.ENFORCE_CAUSALITY
    impulse_response = original_impulse_response;
end
% truncate final samples smaller than 1e-3 of the peak
ir_peak = max(abs(impulse_response));
ir_last  = find(abs(impulse_response)>ir_peak*OP.impulse_response_truncation_threshold, 1, 'last');

voltage = impulse_response(1:ir_last);
t_base = t_base(1:ir_last);

truncation_dB=20*log10(norm(impulse_response(ir_last+1:end))/norm(voltage));
