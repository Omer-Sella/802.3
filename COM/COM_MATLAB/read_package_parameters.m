% end read_sp4_sparam
function param_struct = read_package_parameters(parameter,param_struct)

%With the introduction of Package sections in COM spreadsheet, it makes sense to have a single function that grabs all package parameters from a parameter block
%This block should eventually replace what is in read_ParamConfigFile
%It can be called as:  param = read_package_parameters(parameter, param)

if nargin<2
    %param_struct doesn't need to be passed when building a new package structure
    %it is only needed when appending to regular param structure
    param_struct=struct;
end

param_struct.C_pkg_board = xls_parameter(parameter, 'C_p', true)*1e-9; % C_p in nF (single sided)
param_struct.R_diepad = xls_parameter(parameter, 'R_d', true); % Die source termination resistance  (single sided)

param_struct.a_thru = xls_parameter(parameter, 'A_v', true); % Victim differential peak source output voltage (half of peak to peak)
param_struct.a_fext = xls_parameter(parameter, 'A_fe', true); % FEXT aggressor differential peak source output voltage (half of peak to peak)
param_struct.a_next = xls_parameter(parameter, 'A_ne', true); % NEXT aggressor differential peak source output voltage (half of peak to peak)

param_struct.z_p_tx_cases = xls_parameter(parameter, 'z_p (TX)', true).'; % List of victim transmitter package trace lengths in mm, one per case
[ncases, mele]=size(param_struct.z_p_tx_cases);
if mele ==2
    param_struct.flex=2;
elseif mele==4
    param_struct.flex=4;
elseif mele==1
    param_struct.flex=1;
else
    error('config file syntax error')
end
param_struct.z_p_next_cases = xls_parameter(parameter, 'z_p (NEXT)', true).'; % List of NEXT transmitter package trace lengths in mm, one per case
[ncases1, mele1]=size(param_struct.z_p_next_cases);
if ncases ~= ncases1 || mele ~= mele1
    error('All TX, NEXT, FEXT, Rx cases must agree');
else
end
param_struct.z_p_fext_cases = xls_parameter(parameter, 'z_p (FEXT)', true).'; % List of FEXT transmitter package trace lengths in mm, one per case
[ncases1, mele1]=size(param_struct.z_p_fext_cases);
if ncases ~= ncases1 ||  mele ~= mele1
    error('All TX, NEXT, FEXT, Rx cases must agree');
else
end
param_struct.z_p_rx_cases = xls_parameter(parameter, 'z_p (RX)', true).'; % List of FEXT receiver package trace lengths in mm, one per case
[ncases1, mele1]=size(param_struct.z_p_rx_cases);
if ncases ~= ncases1 ||  mele ~= mele1
    error('All TX, NEXT, FEXT, Rx cases must agree');
else
end
% Table 93A-3 parameters
param_struct.pkg_gamma0_a1_a2 = xls_parameter(parameter, 'package_tl_gamma0_a1_a2', true, [0 1.734e-3 1.455e-4]); %Fitting parameters for package model per unit length. First element is in 1/mm and affects DC loss of package model . Second element is in ns1/2/mm and affects loss proportional to sqrt(f). Third element is in ns/mm and affects loss proportional to f.
param_struct.pkg_tau = xls_parameter(parameter, 'package_tl_tau', true, 6.141e-3); % Package model transmission line delay ns/mm
param_struct.pkg_Z_c = xls_parameter(parameter, 'package_Z_c', true, 78.2).';% Package model transmission line characteristic impedance [ Tx , Rx ]
[ ncases1, mele1]=size(param_struct.pkg_Z_c);%
if   mele ~= mele1
    error('tx rx pairs must have thesame number element entries as TX, NEXT, FEXT, Rx');
else
end
if mele1==2 % fuill in a array if only a 2 element flex package is specified
    for ii=1:ncases
        param_struct.z_p_fext_casesx(ii,:)=  [param_struct.z_p_fext_cases(ii,:)' ;[ 0 ; 0 ]]';
        param_struct.z_p_next_casesx(ii,:)=  [param_struct.z_p_next_cases(ii,:)' ;[ 0 ; 0 ]]';
        param_struct.z_p_tx_casesx(ii,:)=  [param_struct.z_p_tx_cases(ii,:)' ;[ 0 ; 0 ]]';
        param_struct.z_p_rx_casesx(ii,:)=  [param_struct.z_p_rx_cases(ii,:)' ;[ 0 ; 0 ]]';
    end
    param_struct.z_p_fext_cases =  param_struct.z_p_fext_casesx;
    param_struct.z_p_next_cases=  param_struct.z_p_next_casesx;
    param_struct.z_p_tx_cases=  param_struct.z_p_tx_casesx;
    param_struct.z_p_rx_cases=  param_struct.z_p_rx_casesx;
    param_struct.pkg_Z_c=[param_struct.pkg_Z_c' ;[ 100 100 ; 100 100 ]]';
end