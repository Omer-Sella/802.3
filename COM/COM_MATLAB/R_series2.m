function S =R_series2(zref,f,R)
r=ones(1,length(f))*R;
S.Parameters(1,1,:) =  r./(r + 2*zref);
S.Parameters(2,2,:) =  r./(r + 2*zref);
S.Parameters(2,1,:) = (2*zref)./(r + 2*zref);
S.Parameters(1,2,:) =  (2*zref)./(r + 2*zref);
% Sm=sparameters(S.Parameters,f,zref);

function H_tw=Raised_Cosine_Filter(param,f,use_RC)

if use_RC
    H_tw = Tukey_Window(f,param ,param.RC_Start, param.RC_end);% add tw filter;
else
    H_tw=ones(1,length(f));
end