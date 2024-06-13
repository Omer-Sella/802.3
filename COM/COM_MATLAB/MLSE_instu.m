function [MLSE_results] = MLSE_instu(param,alpha,A_s,A_ni,PDF,CDF)
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
    snr_dfe = @(der,PDF,CDF) -10*log10((A_s./CDF_inv_ev(der,PDF,CDF)).^2)+10*log10((L^2-1)/3*qfuncinv(der).^2) ;

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
    %%
    dscale=.05;
    scale=1;
    last_scale_tune=inf;
    scale_tune=inf;
    while abs(scale_tune) >= .1
        istart=-PDF.Min+1;
        scale=scale-dscale;
        PDF_SCALED = scalePDF(PDF,scale);
        cdf_scaled=pdf_to_cdf(PDF_SCALED);
        test_snr=snr_dfe(DER_MLSE_CDF,PDF_SCALED,cdf_scaled.y);
        test_DER=2*(L-1)/L* (CDF_ev(main/(L-1),PDF_SCALED,cdf_scaled.y) );
        scale_tune=(test_DER-DER_MLSE_CDF)/DER_MLSE_CDF;
        if sign(scale_tune) ~= sign(last_scale_tune)
            % scale=scale+dscale % back up
            dscale=-dscale/2;
        end
        last_scale_tune=scale_tune;
    end
    new_com_CDF=10*log10((A_s./CDF_inv_ev(DER0,PDF_SCALED,cdf_scaled.y)).^2);
    delta_com=new_com_CDF-10*log10((A_s./CDF_inv_ev(DER0,PDF,CDF)).^2); 
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
    PDF_SCALED=[];
    cdf_scaled=[];
end

%%
MLSE_results.COM_from_matlab=COM_from_matlab;
MLSE_results.DER_MLSE_Gaussian=DER_MLSE;
MLSE_results.DER_MLSE_CDF=DER_MLSE_CDF;
MLSE_results.sigma_noise=sigma_noise;
MLSE_results.k_DER=k_DER;
MLSE_results.COM_CDF=new_com_CDF;
MLSE_results.delta_com_CDF=delta_com;
MLSE_results.delta_com_Gaussian=delta_com;
MLSE_results.PDF=PDF_SCALED;
MLSE_results.CDF=cdf_scaled.y;
MLSE_results.PDF_scale=scale;

