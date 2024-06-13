function [hisi,tap_coef,hisi_ref]=applyDFEbk(hisi,hisi_ref,idx,tap_bk,curval,bmaxg,dfe_delta)
% hrem=applyDFEbk(hsis,idx,tap_bk,bmaxg)
% applying a bank of DFE at desired location
% hisi: waveform with cursor values;
% idx: starting index of the bank;
% tap_bk: number of taps in bank;
% bmaxg: maxmum coefficient in bank;

if nargin<6
    dfe_delta=0;
end

rng=idx:idx+tap_bk-1;
flt_curval=hisi(rng);

%floor(abs(floatingcursors/sbr(cursor_i))./param.dfe_delta).*param.dfe_delta.*sign(floatingcursors)*sbr(cursor_i);

if dfe_delta~=0
    flt_curval_q=floor(abs(flt_curval/curval)./dfe_delta) .* ...
        dfe_delta.*sign(flt_curval)*curval;
else
    flt_curval_q=hisi(rng);
end

tap_coef=min(abs(flt_curval_q/curval),bmaxg).*sign(flt_curval_q);
hisi(rng)= hisi(rng) - curval*tap_coef;
hisi_ref(rng)=0;

%AJG021820