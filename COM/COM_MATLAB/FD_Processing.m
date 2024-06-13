function [chdata,output_args]=FD_Processing(chdata,output_args,param,OP,SDDp2p,DO_ONCE)
%This function calculates various frequency domain metrics
%Mainly IL_fit, FOM_ILD, ICN, ICN_Fext, and ICN_Next
db = @(x) 20*log10(abs(x));
package_testcase=OP.pkg_len_select(param.package_testcase_i);
if OP.WC_PORTZ
    A_thru = param.a_thru(param.Tx_rd_sel);
    A_fext = param.a_fext(param.Tx_rd_sel);
    A_next = param.a_next(param.Tx_rd_sel);
else
    A_thru = param.a_thru(package_testcase);
    A_fext = param.a_fext(package_testcase);
    A_next = param.a_next(package_testcase);
end
for i=1:param.number_of_s4p_files
    if isequal(chdata(i).type, 'THRU')
        chdata(i).A=A_thru;
        chdata(i).Aicn=A_thru;
    elseif isequal(chdata(i).type, 'FEXT')
        chdata(i).A=A_fext;
        chdata(i).Aicn=param.a_icn_fext;
    elseif isequal(chdata(i).type, 'NEXT')
        chdata(i).A=A_next;
        chdata(i).Aicn=param.a_icn_next;
    end
end
if OP.TDMODE
    for i=1:param.number_of_s4p_files           % freq delta for integration
        chdata(i).delta_f=chdata(i).faxis(11)-chdata(i).faxis(10);
    end
end
if ~DO_ONCE
    return;
end
%Any new output_args fields set in this function should be initialized here as empty
output_args.fitted_IL_dB_at_Fnq = [];
output_args.cable__assembley_loss=[];
output_args.loss_with_PCB=[];
output_args.VIP_to_VMP_IL_dB_at_Fnq=[];
output_args.IL_dB_channel_only_at_Fnq=[];
output_args.VTF_loss_dB_at_Fnq=[];
output_args.IL_db_die_to_die_at_Fnq=[];
output_args.FOM_TDILN=[];
output_args.TD_ILN=[];
output_args.FOM_RILN=[];
output_args.FOM_ILD=[];
%TD_Mode is just a pass through to set the empty values and return
if ~OP.GET_FD
    return;
end
case_number=param.package_testcase_i;
f2=param.f2;
f1=param.f1;
MDFEXT_ICN=0; MDNEXT_ICN=0;
for i=1:param.number_of_s4p_files
    if OP.INCLUDE_FILTER % apply RX filtRaised_Cosine_Filterer
        % Equation 93A-20 %%
        %                     H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(i).faxis./(param.f_r*param.fb));
        f=chdata(i).faxis;
        %
        H_bt=Bessel_Thomson_Filter(param,f,OP.Bessel_Thomson);
        H_bw=Butterworth_Filter(param,f,OP.Butterworth);
        H_RCos=Raised_Cosine_Filter(param,f,OP.Raised_Cosine); % conditionally include the RCos filter for all IR conversion using COM_FD_to_TD
        H_txffe=      Tx_FFE_Filter(param,f,param.Pkg_TXFFE_preset); % RIM 08-18-2022 to add forced TX ffe per package case
        H_r=H_bw.*H_bt.*H_RCos.*H_txffe; % RIM 08-18-2022 to add forced TX ffe per package case
        chdata(i).sdd21=chdata(i).sdd21.*H_r;
        if OP.DISPLAY_WINDOW
            if i==1
                figure(300+param.package_testcase_i);
                subplot(3,1,1)   
                hold on
                plot(chdata(i).faxis/1e9, 20*log10(abs(squeeze(chdata(i).sdd21))), 'k-','linewidth',2, 'Disp',sprintf('VTF (no Tx/Rx eq)')')
                try
                    legend('NumColumns',2)
                    legend('location','south')
                catch
                end
            end
        end
    end
end
for i=1:param.number_of_s4p_files
    if i == 2
        PSXT(1:length(chdata(i).sdd21f))=0;
        MDFEXT(1:length(chdata(i).sdd21f))=0;
        MDNEXT(1:length(chdata(i).sdd21f))=0;
    end
    a=find(chdata(i).faxis(:)>=f2,1,'first');% RIM 01-12-21
    if isempty(a)
        f2=chdata(i).faxis(end);
        index_f2=length(chdata(i).faxis);
    else
        index_f2=a(1);
    end
    b=find(chdata(i).faxis(:)<=f1,1,'last');% RIM 01-12-21
    if isempty(b)
        f1=chdata(i).faxis(1);
        index_f1=1;
    else
        index_f1=b(1);
    end
    % R is the frequency dependent parameter for the sinc function use in the
    % PWF for ICN
    temp_angle=(param.samples_per_ui*param.sample_dt)*pi.*chdata(i).faxis;
    if(chdata(i).faxis(1)==0)
        temp_angle(1)=1e-20;% we don't want to divide by zero
    end
    SINC = sin(temp_angle)./temp_angle;
    PWF_data=SINC.^2;
    PWF_trf=(1+(chdata(i).faxis/chdata(i).ftr).^4).^-1;
    %// bw1=2.613126; bw2=3.4142136; bw3=2.613126;
    fr=param.f_r*param.fb;
    PWF_rx=(1+(chdata(i).faxis/fr).^8).^-1;
    PWF_highpass=1;
    % Equation 93A-57 %
    PWF=PWF_data.*PWF_trf.*PWF_rx.*PWF_highpass; % power weight function
    % freq delta for integration
    chdata(i).delta_f=chdata(i).faxis(11)-chdata(i).faxis(10);
    % from ba spec, this is basically ICN
    faxis_GHz = chdata(i).faxis/1e9;
    if isequal(chdata(i).type, 'THRU')
                [ILD_magft chdata(i).fit_orig] = get_ILN(chdata(i).sdd21f(index_f1:index_f2), chdata(i).faxis(index_f1:index_f2));
        % find fitted loss values by interpolation using full data, no indexing - Adee 2022-08-28
        [~, chdata(i).fit_orig] = get_ILN(chdata(i).sdd21f, chdata(i).faxis);
        fit_loss = interp1(chdata(i).faxis, -chdata(i).fit_orig, 1/param.ui/2);
        chdata(i).fit_ILatNq = fit_loss;
        output_args.fitted_IL_dB_at_Fnq = fit_loss;
        IL_interp = interp1(chdata(i).faxis, -20*log10(abs(chdata(i).sdd21f)), 1/param.ui/2);
        chdata(i).ILatNq = IL_interp;
        if OP.include_pcb
            cable_loss = interp1(chdata(i).faxis, -20*log10(abs(chdata(i).sdd21_orig)), 1/param.ui/2);
            loss_with_PCB = interp1(chdata(i).faxis, -20*log10(abs(chdata(i).sdd21_raw)), 1/param.ui/2);
            output_args.cable__assembley_loss=cable_loss;
            output_args.loss_with_PCB=loss_with_PCB;
        end
        Nq_loss=chdata(i).ILatNq;
        output_args.IL_dB_channel_only_at_Fnq=Nq_loss;
        % time domain ref RR = complex fit pulse
        if OP.COMPUTE_TDILN || OP.COMPUTE_RILN
            [ILD chdata(i).fit, TD_ILN ] = get_ILN_cmp_td(chdata(i).sdd21f(index_f1:index_f2), chdata(i).faxis(index_f1:index_f2),OP,param,chdata(i).A);
            FOM_TDILN = TD_ILN.SNR_ISI_FOM_PDF;
            FOM_ILN_complex= TD_ILN.FOM;
        end
        if OP.COMPUTE_TDILN || OP.COMPUTE_RILN
            [ILD chdata(i).fit, TD_ILN ] = get_ILN_cmp_td(chdata(i).sdd21f(index_f1:index_f2), chdata(i).faxis(index_f1:index_f2),OP,param,chdata(i).A);
            FOM_TDILN = TD_ILN.SNR_ISI_FOM_PDF;
            FOM_ILN_complex= TD_ILN.FOM;
        end
        if OP.COMPUTE_TDILN
            output_args.FOM_TDILN=FOM_TDILN;
            output_args.TD_ILN=TD_ILN; % struct
        end
        if OP.COMPUTE_RILN
            % Get RIL, RILN, and TD_RILN
            [RIL_struct]= capture_RIL_RILN(chdata);
            FOM_RILN=sqrt(chdata(i).delta_f/(param.f2-param.f1)*sum( PWF(index_f1:index_f2-1).*RIL_struct.RILN_dB(index_f1:index_f2-1)'.^2));
            output_args.FOM_RILN=FOM_RILN;           
            %---start. plotting ILN based on ILD and RILN % Hansel 10/18/2021
            plot_tdomain_debug= 0; % must have OP.COMPUTE_TDILN = 1 to use
            if plot_tdomain_debug== 1
                figure(988); set(gcf,'Tag','COM')
                ax_1= subplot(3,1,1);
                plot(TD_ILN.REF.t*1e9, TD_ILN.REF.PR*1e3,'disp','ref');
                hold on;
                plot(TD_ILN.FIT.t*1e9, TD_ILN.FIT.PR*1e3,'disp','fit');
                hold on;
                plot(TD_ILN.t*1e9, TD_ILN.ILN*1e3, 'k','disp','ref - fit (IL noise)');
                ylim([min(TD_RILN.ILN) max(TD_RILN.ILN)]*1e3);
                grid on;
                box on;
                legend('REF', 'FIT', 'TD\_ILN: ref - fit (IL noise)');
                xlabel('Time [nsec]');
                ylabel('Pulse Response [mV]');
                
                ax_2= subplot(3,1,2);
                plot(TD_RILN.REF.t*1e9, TD_RILN.REF.PR*1e3);
                hold on;
                plot(TD_RILN.t*1e9, TD_RILN.ILN*1e3, 'r');
                ylim([min(TD_RILN.ILN) max(TD_RILN.ILN)]*1e3);
                grid on;
                box on;
                legend('REF', 'TD\_RILN');
                xlabel('Time [nsec]');
                ylabel('Pulse Response [mV]');
                ax_3= subplot(3,1,3);
                plot(TD_ILN.t*1e9, TD_ILN.ILN*1e3, 'k');
                hold on;
                plot(TD_RILN.t*1e9, TD_RILN.ILN*1e3, 'r');
                ylim([min(TD_RILN.ILN) max(TD_RILN.ILN)]*1e3);
                grid on;
                box on;
                legend( 'TD\_ILN: ref - fit (IL noise)', 'TD\_RILN');
                xlabel('Time [nsec]');
                ylabel('Pulse Response [mV]');
                
                linkaxes([ax_1, ax_2, ax_3], 'x');
                ax_1.XLim = [0 max(TD_RILN.t)*1e9 ];
            end
            %---end. plotting ILN based on ILD and RILN
        end
        % Equation 93A-56 %
        FOM_ILD=sqrt(chdata(i).delta_f/(param.f2-param.f1)*sum( PWF(index_f1:index_f2).*ILD_magft.^2));
        output_args.FOM_ILD=FOM_ILD;
        if OP.DEBUG
            if OP.DISPLAY_WINDOW                
                figure(300+case_number);
                set(gcf,'Tag','COM')
                screen_size=get(0,'ScreenSize');
                pos = get(gcf, 'OuterPosition');
                set(gcf, 'Name', [sprintf('%.3gdB IL Channel: ',Nq_loss) 'Raw frequency-domain data'], 'OuterPosition', ...
                    screen_size([3 4 3 4]).*[0 1 0 0] + pos([3 4 3 4]).*[0 -2 1 2] ...IL fit
                    - (case_number-1)*[0 20 0 0]);
                subplot(3,1,1)
                title('Losses')
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd21f))), 'b', 'LineWidth', 3, 'Disp','IL passed s-params')
                hold on
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd21p_nodie))), 'm-', 'Disp','IL die to die (w pkg/brds)')
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd21p))), 'b-', 'Disp','IL dB VIP to VMP')               
                ylim(get(gca, 'ylim'));
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd11))),'c','Disp','RL11')
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd22))),'m','Disp','RL22')
                subplot(3,1,3)
                plot(faxis_GHz(index_f1:index_f2), ILD_magft,'Disp','ILD')
                if OP.PLOT_CM
                    if case_number ==1
                        h350=figure(350);set(gcf,'Tag','COM')
                        screen_size=get(0,'ScreenSize');
                        pos = get(gcf, 'OuterPosition');
                        set(gcf, 'OuterPosition',screen_size([3 4 3 4]).*[1  1 0 0] + pos([3 4 3 4]).*[-2 -2 2 1])
                        movegui(gcf,'center');
                        htabgroup350 = uitabgroup(h350);
                        htab1 = uitab(htabgroup350, 'Title', 'CM Through Losses');
                        hax1 = axes('Parent', htab1);
                        set(h350,'CurrentAxes',hax1)
                        hold on
                        set(gcf,'Tag','COM')
                        screen_size=get(0,'ScreenSize');
                        pos = get(gcf, 'OuterPosition');
                        title('IL & CM Losses')
                        base=strrep(chdata(i).base,'_',' ');
                        plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd21f))), 'b', 'LineWidth', 3, 'Disp',  [ 'sdd21(IL) TP0-TP5 ' base])
                        plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdc21_raw))),   'LineWidth', 2, 'Disp',[ 'sdc21 TP0-TP5 '  base])
                        ylabel('dB')
                        xlabel('GHz')
                        legend show
                        legend('Location','eastoutside')
                        hold on
                        grid on
                        if param.number_of_s4p_files > 1
                            htab2 = uitab(htabgroup350, 'Title', 'CM Crosstalk Losses');
                            hax2 = axes('Parent', htab2);
                            htab3 = uitab(htabgroup350, 'Title', 'Crosstalk Losses');
                            hax3 = axes('Parent', htab3);
                        end
                        
                    end
                end
            else
                display(['Insertion Loss at Nyquist = ', num2str(chdata(i).ILatNq)])
            end
        end
    else % NEXT or FEXT
        if isequal(chdata(i).type, 'FEXT')
            MDFEXT=sqrt(abs(chdata(i).sdd21f).^2+MDFEXT.^2); % power sum xtk
            MDFEXT_ICN=sqrt(2*chdata(i).delta_f/param.f2*sum( chdata(i).Aicn^2*PWF(index_f1:index_f2).*abs(MDFEXT(index_f1:index_f2)).^2)); %eq 46
            output_args.MDFEXT_ICN_92_47_mV=MDFEXT_ICN*1000;
        elseif isequal(chdata(i).type, 'NEXT')
            MDNEXT=sqrt(abs(chdata(i).sdd21f).^2+MDNEXT.^2); % power sum xtk
            MDNEXT_ICN=sqrt(2*chdata(i).delta_f/param.f2*sum( chdata(i).Aicn^2*PWF(index_f1:index_f2).*abs(MDNEXT(index_f1:index_f2)).^2)); %eq 47
            output_args.MDNEXT_ICN_92_46_mV=MDNEXT_ICN*1000;
        end
        PSXT=sqrt((abs(chdata(i).sdd21f)*chdata(i).Aicn).^2+PSXT.^2); % power sum xtk
        ICN=sqrt(2*chdata(i).delta_f/param.f2*sum( PWF(index_f1:index_f2).*abs(PSXT(index_f1:index_f2)).^2));
        output_args.ICN_mV=ICN*1000;
        ICN_test=norm([MDFEXT_ICN MDNEXT_ICN]);
        if  OP.PLOT_CM && OP.DISPLAY_WINDOW
            if case_number ==1
                %                             htab2 = uitab(htabgroup350, 'Title', 'CM Crosstalk Losses');
                %                             hax2 = axes('Parent', htab2);
                set(h350,'CurrentAxes',hax2)
                hold on
                title('CM Losses')
                base=strrep(chdata(i).base,'_',' ');
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdc21_raw))), 'LineWidth', 2, 'Disp',[ 'sdc21 TP0-TP5 ' base ])
                legend('Location','eastoutside')
                hold on
                grid on
                set(h350,'CurrentAxes',hax3)
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd21_raw))), 'LineWidth', 2, 'Disp',[ 'sdd21 TP0-TP5 ' base ])
                legend('Location','eastoutside')
                hold on
                grid on
            end
        end
    end
end  % for loop
ICN_test=norm([MDFEXT_ICN MDNEXT_ICN]);
if OP.DEBUG && OP.DISPLAY_WINDOW
    figure(300+case_number);set(gcf,'Tag','COM');
    if param.number_of_s4p_files > 1
        scale=1/chdata(2).Aicn; %chdata(i).sdd21f not scalled
        subplot(3,1,1)
        hold on
        plot(faxis_GHz, 20*log10(abs(PSXT*scale)),'r','Disp','PSXTK')
        icrxi=find(chdata(i).faxis >=param.fb/2,1,'first');
        subplot(3,1,2)
        grid on
        ILtemp=20*log10(abs(chdata(1).sdd21f));
        IL4ICR=interp1(chdata(1).faxis,ILtemp,chdata(i).faxis);
        scale=1/chdata(2).Aicn; % chdata(i).sdd21f not scalled
        ICR=-20*log10(abs(PSXT*scale))+IL4ICR;
        semilogx(faxis_GHz, ICR,'Disp', 'ICR')
        hold on
        stem(faxis_GHz(icrxi), ICR(icrxi),'g', 'disp', 'f_{Baud}/2')
    end
    subplot(3,1,1)
    title([param.base ' Losses']); ylabel('dB'); xlabel('GHz')
    grid on;  legend show
    subplot(3,1,2)
    title([param.base ' ICR']); ylabel('dB'); xlabel('GHz')
    ylim([0 80])
    xlim([.1 100])
    grid on; %legend show
    subplot(3,1,3)
    title([param.base ' ILD']); ylabel('dB'); xlabel('GHz')
    ylim([-3 3])
    grid on; legend show
end
% find loss values by interpolation using full data, no indexing - Adee 2022-08-28
total_loss = interp1(chdata(i).faxis, -20*log10(abs(chdata(1).sdd21)), 1/param.ui/2);
d2d_loss = interp1(chdata(i).faxis, -20*log10(abs(chdata(1).sdd21p_nodie)), 1/param.ui/2); 
VIP_VMP_loss = interp1(chdata(i).faxis, -20*log10(abs(chdata(1).sdd21p)), 1/param.ui/2); 
output_args.VTF_loss_dB_at_Fnq=total_loss;
output_args.IL_db_die_to_die_at_Fnq=d2d_loss;
output_args.VIP_to_VMP_IL_dB_at_Fnq=VIP_VMP_loss;