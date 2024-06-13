function [tap_loc,tap_coef,hisi,b]=floatingDFE( hisi, N_b, N_bf, N_bg, N_bmax, bmaxg, curval, dfe_delta)

% hisi = postcursor isi
% N_b = number of fixed dfe taps (before floating taps begin)
% N_bf = number of floating taps per group
% N_bg = nubmber of groups
% N_bmax = max tap number that can be used for floating tap
% bmaxg = max tap strength for floating taps
% curval = value of the cursor


if nargin<8, dfe_delta=0;end


tap_coef=zeros(1,length(hisi));
b=zeros(1,length(hisi));


[tap_loc]=findbankloc(hisi,N_b+1,N_bmax,N_bf,curval,bmaxg,N_bg);

%Apply DFE to all taps
flt_curval=hisi(tap_loc);
if dfe_delta~=0
    flt_curval_q=floor(abs(flt_curval/curval)./dfe_delta) .* ...
        dfe_delta.*sign(flt_curval)*curval;
else
    flt_curval_q=hisi(tap_loc);
end
applied_coef=min(abs(flt_curval_q/curval),bmaxg).*sign(flt_curval_q);
hisi(tap_loc)= hisi(tap_loc) - curval*applied_coef;
tap_coef(tap_loc)=applied_coef;



tap_loc=sort(tap_loc,'ascend');
b(tap_loc)=bmaxg;