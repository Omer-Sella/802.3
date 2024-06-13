function [ FIR t] =get_RAW_FIR(H,f,OP,param)
H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*f./(0.75*param.fb));
if ~iscolumn(H), H=H.';end
if ~iscolumn(H_r), H_r=H_r.';end
H=H(:).*H_r;
[FIR, t, ~,~] = s21_to_impulse_DC(H ,f, param.sample_dt, OP) ;
% SBR=filter(ones(1, param.samples_per_ui), 1, FIR);