function [ sigma_NE , sigma_HP] = get_sigma_noise( H_ctf, param, chdata, sigma_bn )
% for Rx calibratrion only
% the FEXT channel for calibration basically a DC connection unlike normal
% FEXT channels which are nearly open at DC channels
H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(2).faxis/(param.f_r*param.fb));
idxfbby2=find( chdata(2).faxis(:) >= param.fb/2, 1);
if size(chdata,2) >= 2
    Hnoise_channel=chdata(2).sdd21;% rx package is already included tx is not
else
    Hnoise_channel=1;
end
f=chdata(2).faxis;
f_hp=param.f_hp;
if f_hp ~=0 % param.f_hp is a key indicating that Tx bbn is used as in clause 162 
    H_hp=(-1j*f./f_hp)./(1+1j*f./f_hp);
else
    H_hp=ones(1,length(f));
end
%% Equation 93A-47 or 162-12 
H_np=Hnoise_channel.*H_ctf.*H_r.*H_hp;

%% Equation 93A-48 or 162-14%%
sigma_NE = sigma_bn*sqrt(mean(abs(H_np(1:idxfbby2).^2)));
sigma_HP = sigma_bn*(mean(abs(H_hp(1:idxfbby2).^2)));