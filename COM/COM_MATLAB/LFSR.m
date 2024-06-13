function[seq c]=LFSR(s,t)
%s=initial state of LFSR, you can choose any lenght of LFSR
%Instruction:========== 
%Save LFSR.m in your current directory and type following
%on Command window for simulating 5 bit LFSR with tap [5 2]
%---------------------
%>>s=[1 1 0 0 1] 
%>>t=[5 2]
%>>[seq c] =LFSR(s,t) 
%---------------------------
%seq = generated sequence
%c will be matrix containing the states of LFSR raw wise
% 
%-----------------------------------------------------------
%If any doubt, confusion or feedback please contact me
% NIKESH BAJAJ 
% bajaj.nikkey@gmail.com (+91-9915522564)
% Asst. Professor at Lovely Profesional University
% Masters from Aligarh Muslim University,INDIA 
%--------------------------------------------------
n=length(s);
c(1,:)=s;
m=length(t);
for k=1:2^n-2;
b(1)=xor(s(t(1)), s(t(2)));
if m>2;
    for i=1:m-2;
    b(i+1)=xor(s(t(i+2)), b(i));
    end
end
j=1:n-1;
s(n+1-j)=s(n-j);
s(1)=b(m-1);
c(k+1,:)=s;
end
seq=c(:,n)';

function [ dataout ] = pam( data )
% mapping data usng Grey Coding
for i=1:2:floor(length(data)/2)*2
    if data(i:i+1)==[ -1 -1 ]
        dataout(ceil(i/2)) = -1;
    elseif data(i:i+1)==[ -1 1 ]
        dataout(ceil(i/2)) = -1/3;
    elseif data(i:i+1)==[ 1 1 ]
        dataout(ceil(i/2)) = 1/3;
    elseif data(i:i+1)==[ 1 -1 ]
        dataout(ceil(i/2)) = 1;
    end
end
function RILN_TD_struct=RILN_TD(sdd21,RIL,faxis_f2,OP,param,A_T)
db = @(x) 20*log10(abs(x));
disp('computing TD_RILN...')
sdd21=squeeze(sdd21);
if  iscolumn(sdd21)
    sdd21=sdd21.';
end
RIL=squeeze(RIL);
if  iscolumn(RIL)
    RIL=RIL.';
end
print_for_codereview=1;
if exist('OP','var')
    X=sinc(faxis_f2*param.ui)*param.ui*1e9;
    
    H_bt=Bessel_Thomson_Filter(param,faxis_f2,1);
    H_bw=Butterworth_Filter(param,faxis_f2,1);
    H_t = exp(-(pi*faxis_f2/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
    H_tw=Tukey_Window(faxis_f2,param);
    H_tw=ones(1,length(faxis_f2) );
    [RILN_TD_struct.REF.FIR, ...
        RILN_TD_struct.REF.t, ...
        RILN_TD_struct.REF.causality_correction_dB, ...
        RILN_TD_struct.REF.truncation_dB] = s21_to_impulse_DC(sdd21.*H_bw.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    RILN_TD_struct.REF.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.REF.FIR);
    [RILN_TD_struct.FIT.FIR, ...
        RILN_TD_struct.FIT.t, ...
        RILN_TD_struct.FIT.causality_correction_dB, ...
        RILN_TD_struct.FIT.truncation_dB] = s21_to_impulse_DC(RIL.*H_bw.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    RILN_TD_struct.FIT.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.FIT.FIR);   
    ipeak=find(RILN_TD_struct.REF.PR==max(RILN_TD_struct.REF.PR),1,'first');
    NrangeUI=1000;
    range_end=min(ipeak+param.samples_per_ui*NrangeUI,min(length(RILN_TD_struct.FIT.FIR), length(RILN_TD_struct.REF.FIR) ) ) ; 
    range=ipeak:range_end;
    RILN_TD_struct.ILN=RILN_TD_struct.FIT.PR(range)-RILN_TD_struct.REF.PR(range);
    RILN_TD_struct.t=RILN_TD_struct.FIT.t(range);
    RILN_TD_struct.FOM=-inf;
    RILN_TD_struct.FOM_PDF=-inf;
    rms_fom=-inf;
    for im=1:param.samples_per_ui
        RILN_TD_struct.FOM=max(RILN_TD_struct.FOM, norm( RILN_TD_struct.ILN(im:param.samples_per_ui:end)));
        [ pdf ] = get_pdf_from_sampled_signal(  RILN_TD_struct.ILN(im:param.samples_per_ui:end), param.levels, OP.BinSize ,0);
        rms=sqrt(pdf.y*pdf.x(:).^2)*sqrt(2);
        cdf=pdf; cdf.y=cumsum(pdf.y);
        %         cursors = d_cpdf(OP.BinSize,param.a_thru*[-1 1], [1 1]/2);
        %         signal_and_isi_pdf = conv_fct(cursors, pdf);
        %         cdf=signal_and_isi_pdf; cdf.y=cumsum(signal_and_isi_pdf.y);
        if print_for_codereview % remove once all checked out
            h=figure(191);set(gcf,'Tag','COM');
            semilogy(-cdf.x,cdf.y);
%             xlim ([0,-cdf.x(1)])
            ylim([param.specBER 1]);title ('CDF of RILN')
            hold on
        end
        if rms>rms_fom
            rms_fom=rms;
            RILN_TD_struct.FOM_PDF= -cdf.x(find(cdf.y >= param.specBER, 1, 'first'));
            RILN_TD_struct.PDF=pdf;
        end
    end
    pdf_from_norm=normal_dist(RILN_TD_struct.FOM, 7 , OP.BinSize);
    RILN_TD_struct.SNR_ISI_FOM=db(RILN_TD_struct.FIT.PR(ipeak)/RILN_TD_struct.FOM);
    RILN_TD_struct.SNR_ISI_FOM_PDF=db(RILN_TD_struct.FIT.PR(ipeak)/RILN_TD_struct.FOM_PDF);
    if print_for_codereview % remove once all checked out
        figure(9003);set(gcf,'Tag','COM');
        plot(RILN_TD_struct.t,RILN_TD_struct.ILN,'disp','td riln')
        hold on
        plot(RILN_TD_struct.FIT.t,RILN_TD_struct.FIT.PR,'disp','fit')
        plot(RILN_TD_struct.REF.t,RILN_TD_struct.REF.PR,'disp','ref')
        hold off
        fprintf('SNR ISI FOM rms = %g dB;   SNR ISI FOM PDF = %g dB\n',RILN_TD_struct.SNR_ISI_FOM,RILN_TD_struct.SNR_ISI_FOM_PDF)
        figure(9004);set(gcf,'Tag','COM');
        semilogy(RILN_TD_struct.PDF.x,RILN_TD_struct.PDF.y,'disp','actual PDF')
        hold on
        semilogy(pdf_from_norm.x,pdf_from_norm.y,'disp','PDF using Gaussian assumed PDF');
        ylim([param.specBER max([RILN_TD_struct.PDF.y pdf_from_norm.y])]);title ('Compare actual PDF to Gaussian')
        grid on
        legend('show')
    end
end