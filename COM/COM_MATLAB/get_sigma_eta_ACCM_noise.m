function [sigma_N] = get_sigma_eta_ACCM_noise(chdata,param,H_sy,H_r,H_ctf)
% 11-25-2020 correct integratation limits, should be infinity or highest specified frequency
% H_sy - PSD for Tx power delivery, not normally used and set to a vector 1's
% H_r - receiver filter, Butterworth
% H_ctf - total gain of CTLE and low freq filtering
% H_dc - the common mode channel gain
% param.eta_0 -input referred Rx noise
% param.AC_CM_RMS - AC CM source before Tx series source resistor.
% param.ACCM_MAX_Freq - Max frequency for ACCM source intergration
%% Equation 93A-35 - independent of FFE setting %%
sigma_N1 = sqrt(param.eta_0*sum( abs(H_sy(2:end) .* H_r(2:end) .* H_ctf(2:end) ).^2 .* diff(chdata(1).faxis)/1e9));% eta_0 is V^2/Ghz i.e. /1
if sum(param.AC_CM_RMS) ~= 0
     sigma_ACCM=0;
    f_int= chdata(1).faxis( chdata(1).faxis<=param.ACCM_MAX_Freq );
    for i=1:length(chdata)       
        H_dc=abs(squeeze(chdata(i).sdc21));
        sigma_ACCM_acc= sqrt( 2*param.AC_CM_RMS_TX^2.*sum( abs(H_sy(2:length(f_int)) .* H_r(2:length(f_int)) .* H_ctf(2:length(f_int)).*H_dc(2:length(f_int)) ).^2 .* diff(f_int))/f_int(end) );
        sigma_ACCM=norm([sigma_ACCM_acc,sigma_ACCM]);
    end
    sigma_N=norm([sigma_N1,sigma_ACCM]);
else
    sigma_N=sigma_N1;
end
%%