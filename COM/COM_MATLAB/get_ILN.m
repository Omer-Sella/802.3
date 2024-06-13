function [ILN, efit]= get_ILN(sdd21,faxis_f2)
% used for FD IL fitting
% sdd21 us a complex insertion loss
db = @(x) 20*log10(abs(x));
sdd21=squeeze(sdd21);
if  iscolumn(sdd21)
    sdd21=sdd21.';
end
fmbg=[ones(length(faxis_f2),1).*transpose(abs(sdd21))  transpose(sqrt(faxis_f2)).*transpose(abs(sdd21))  transpose(faxis_f2).*transpose(abs(sdd21)) transpose(faxis_f2.^2).*transpose(abs(sdd21)) ];
warning('off','MATLAB:nearlySingularMatrix');
LGw=transpose(abs(sdd21).*db(sdd21));
alpha = ((fmbg'*fmbg)^-1)*fmbg'*LGw;
efit=(alpha(1)+alpha(2).*sqrt(faxis_f2)+alpha(3).*faxis_f2 +faxis_f2.^2.*alpha(4)   );
ILN = db(sdd21)-efit;