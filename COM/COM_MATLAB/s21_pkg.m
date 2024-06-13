
function [s21p,SCH,sigma_ACCM_at_tp0] = s21_pkg(chdata, param, OP, channel_number,mode,include_die)
% concatenates package reflections with s21,s11,and s22 with spec return loss (gammas)
% faxis is the frequency array
% s21, s11, s22 are the corresponding array of differential parameters
% s21p includes the VFT and Tx filter if include_die=1
if nargin<6
    include_die=1;
end
if nargin<5
    mode='dd';
end

s21=chdata.(['s' mode '21_raw']);
s12=chdata.(['s' mode '12_raw']);
s11=chdata.(['s' mode '11_raw']);
s22=chdata.(['s' mode '22_raw']);
faxis=chdata.faxis;
channel_type=chdata.type;

if strcmpi(mode,'dd')
    s11=s11*param.kappa1;
    s22=s22*param.kappa2;
end


Z0=param.Z0;
%sigma_ACCM_at_tp0 is only used when mode=DC
sigma_ACCM_at_tp0=0;

% The following three parameters have possibly different valuesF for TX and
% RX (so can be 2-element vectors).
R_diepad = param.R_diepad;
