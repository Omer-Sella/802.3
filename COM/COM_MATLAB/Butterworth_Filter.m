function H_bw=Butterworth_Filter(param,f,use_BW)

if use_BW
    H_bw = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*f./(param.fb_BW_cutoff*param.fb));
else
    H_bw=ones(1,length(f));
end
function [CDF_ev] = CDF_ev(val,PDF,CDF)
index=find(PDF.x >= -val,1,'first');
CDF_ev=CDF(index);
function [CDF_inv_ev] = CDF_inv_ev(val,PDF,CDF)
index=find(CDF >= val,1,'first');
CDF_inv_ev=PDF.x(index);