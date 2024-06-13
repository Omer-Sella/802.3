function pdf=Init_PDF_Fast( EmptyPDF, values, probs)
%  p=cpdf(type, ...)
%
% CPDF is a probability mass function for discrete distributions or an
% approxmation of a PDF for continuous distributions.
%
% cpdf is internally normalized so that the sum of probabilities is 1
% (regardless of bin size).

% Internal fields:
% Min: *bin number* of minimum value.
% BinSize: size of PDF bins. Bin center is the representative value.
% Vec: vector of probabilities per bin.

pdf=EmptyPDF;

rounded_values_div_binsize=round(values/pdf.BinSize);
%values=pdf.BinSize*rounded_values_div_binsize;

% %speed up for small values round to 0 (because they are all much smaller than binsize)
% if all(values==0)
%     return;
% end
% 
% %speed up for all values rounded to the same bin
% %The output pdf is the same as the
% %empty pdf, but the x value is non-zero (but still scalar)
% if all(values==values(1))
%     pdf.Min=rounded_values_div_binsize(1);
%     pdf.x=values(1);
%     return;
% end
% 
% %The code below requires that values is
% %sorted.  Generally this should be true, but check to be sure
% if ~issorted(values)
%     [values,si]=sort(values);
%     rounded_values_div_binsize=rounded_values_div_binsize(si);
%     probs=probs(si);
% end


%pdf.x=values(1):pdf.BinSize:values(end);
pdf.x=pdf.BinSize*rounded_values_div_binsize(1):pdf.BinSize:pdf.BinSize*rounded_values_div_binsize(end);
pdf.Min=rounded_values_div_binsize(1);

pdf.y=zeros(size(pdf.x));
%The rounded values divided by binsize will reveal the bin number if
%pdf.Min is subtracted from it
bin_placement=rounded_values_div_binsize-pdf.Min+1;
%Can avoid one addition by inserting the first probability
%actually helps when calling this 2 million times
pdf.y(bin_placement(1))=probs(1);
for k=2:length(values)
    pdf.y(bin_placement(k)) = pdf.y(bin_placement(k))+probs(k);
end


%Have already ensured that sum(pdf.y)=1
%pdf.y=pdf.y/sum(pdf.y);

% if any(~isreal(pdf.y)) || any(pdf.y<0)
%     error('PDF must be real and nonnegative');
% end

% pMax=pdf.Min+length(pdf.y)-1;
% pdf.x = values(1):pdf.BinSize:pMax*pdf.BinSize;
function [MLSE_results] = MLSE(param,alpha,A_s,A_ni,PDF,CDF)
% OP.MLSE= 1 ... COM and VEC will be adjusted with MLSE CDF  
% OP.MLSE= 2 ... COM and VEC will be adjusted with MLSE Gaussian assumptions  
% Based on oif2022.580.00 / IEEE802. shakiba_3dj_01_230116 by Hossein Shakiba

qfuncinv = @(x) sqrt(2)*erfcinv(2*x);
qfunc = @(x) 0.5*erfc(x/sqrt(2));

%% step 0
    COM_from_matlab=20*log10(A_s/A_ni);
    L=param.levels;
    DER0=param.specBER;
%% step 1 from slide 6/5
    A_peak=(L-1)*A_s; % slide 6 A_s is main in appendix a
    main=A_peak;
    k_DER=qfuncinv(param.specBER);
    sigma_noise=sqrt(sum(PDF.y.*PDF.x.^2));
    SNR_dB=10*log10( 1/3*(L+1)/(L-1)*(A_peak^2)/sigma_noise^2) ;
    COM=SNR_dB-10*log10((L^2-1)/3*k_DER^2);
%     sprintf('COM from Matlab %g dB\n COM from slide 6 using Gaussian asumptions %g dB\n', COM_from_matlab ,COM)
if A_s >= A_ni
%% step 2 slide 10/8
    SNR_DFE=1/3*(L+1)/(L-1)*(A_peak^2)/sigma_noise^2;
%% step 2 slide 10/8
%     DER_DFE=    2/ ( L/(L-1) -qfunc(   (1-2*alpha)*main/(L-1)/sigma_noise )   )*(qfunc(main/(L-1)/sigma_noise));
%     DER_DFE_CDF=2/ ( L/(L-1)-CDF_ev(   (1-2*alpha)*main/(L-1),PDF,CDF )       )*CDF_ev((main/(L-1)),PDF,CDF);
%% step 3 side 11/9
    j=1:200;
    DER_MLSE=2*sum( j .* ((L-1)/L).^j .* qfunc( sqrt(1+(j-1)*(1-alpha)^2+alpha^2).* main/((L-1)*sigma_noise )       ));
    DER_MLSE_CDF=0; jj=1;
    DER_delta = inf;
    while DER_delta > .001
        last_DER_MLSE_CDF=DER_MLSE_CDF;
        DER_MLSE_CDF=2*( jj .* ((L-1)/L).^jj .* CDF_ev( sqrt(1+(jj-1)*(1-alpha)^2+alpha^2).* main/((L-1) ),PDF,CDF       ))+DER_MLSE_CDF;
        DER_delta= 1-last_DER_MLSE_CDF/DER_MLSE_CDF;
        jj=jj+1;
    end
%% step 4 slide 12/10
    SNR_DFE_eqivalent=SNR_DFE*(...
        (L-1)*sigma_noise/main *    qfuncinv(...
        1/2 *DER_MLSE*(L/(L-1) - qfunc((1-2*alpha)*main/(L-1)*sigma_noise )) ...
        )   ...
        )^2;
    SNR_DFE_eqivalent_CDF=SNR_DFE*(...
        (L-1)/main *    CDF_inv_ev(...
        1/2 *DER_MLSE_CDF*(L/(L-1) - CDF_ev((1-2*alpha)*main/(L-1),PDF,CDF )) ...
        ,PDF, CDF )   ...
        )^2;

%% step 5 slide 13/11
delta_com=10*log10(SNR_DFE_eqivalent/SNR_DFE);
delta_com_CDF=10*log10(SNR_DFE_eqivalent_CDF/SNR_DFE);
new_com_CDF=COM_from_matlab+delta_com_CDF;
else
    warning('MLSE not applied because there is more noise than signal')
    DER_MLSE=[];
    DER_MLSE_CDF=[];
    SNR_DFE_eqivalent=[];
    SNR_DFE_eqivalent_CDF=[];
    new_com_CDF=COM_from_matlab;
    delta_com_CDF=0;
    delta_com=0;
    SNR_DFE=[];
end

%%
MLSE_results.COM_from_matlab=COM_from_matlab;
MLSE_results.SNR_DFE=SNR_DFE;
MLSE_results.DER_MLSE_Gaussian=DER_MLSE;   
MLSE_results.DER_MLSE_CDF=DER_MLSE_CDF;
MLSE_results.sigma_noise=sigma_noise;
MLSE_results.SNR_dB=SNR_dB ;
MLSE_results.SNR_DFE_eqivalent_Gaussian=SNR_DFE_eqivalent;
MLSE_results.SNR_DFE_eqivalent_CDF=SNR_DFE_eqivalent_CDF;
MLSE_results.COM_Gaussian=new_com_CDF;
MLSE_results.COM_CDF=new_com_CDF;
MLSE_results.k_DER=k_DER;
MLSE_results.delta_com_CDF=delta_com_CDF;
MLSE_results.delta_com_Gaussian=delta_com;
