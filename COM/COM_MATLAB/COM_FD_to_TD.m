function chdata=COM_FD_to_TD(chdata,param,OP)

% get impulse responses which in interim step between equation for X(f) and
% H^(k)(t) without TX FFE or CTLE. These will we added later.
case_number=param.package_testcase_i;
for i=1:param.number_of_s4p_files
    %  RIM 2-01-2023 moved to FD_Processing
%     if OP.INCLUDE_FILTER % apply RX filtRaised_Cosine_Filterer
%         % Equation 93A-20 %%
%         %                     H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(i).faxis./(param.f_r*param.fb));
%         f=chdata(i).faxis;
%         %
%         H_bt=Bessel_Thomson_Filter(param,f,OP.Bessel_Thomson);
%         H_bw=Butterworth_Filter(param,f,OP.Butterworth);
%         H_RCos=Raised_Cosine_Filter(param,f,OP.Raised_Cosine); % conditionally include the RCos filter for all IR conversion using COM_FD_to_TD
%         H_txffe=      Tx_FFE_Filter(param,f,param.Pkg_TXFFE_preset); % RIM 08-18-2022 to add forced TX ffe per package case
%         H_r=H_bw.*H_bt.*H_RCos.*H_txffe; % RIM 08-18-2022 to add forced TX ffe per package case
%         chdata(i).sdd21=chdata(i).sdd21.*H_r;
%         if OP.DISPLAY_WINDOW
%             if i==1
%                 figure(300+param.package_testcase_i);
%                 subplot(3,1,1)
%                 plot(chdata(i).faxis/1e9, 20*log10(abs(squeeze(chdata(i).sdd21))), 'k-','linewidth',2, 'Disp','VTF (no Tx/Rx eq)')
%                 try
%                     legend('NumColumns',2)
%                     legend('location','south')
%                 catch
%                 end
%             end
%         end
%     end
    [chdata(i).uneq_imp_response, ...
        chdata(i).t, ...
        chdata(i).causality_correction_dB, ...
        chdata(i).truncation_dB] = s21_to_impulse_DC(chdata(i).sdd21 ,chdata(i).faxis, param.sample_dt, OP) ;
    if ~OP.RX_CALIBRATION || i==1 % DC (common to differentail model is not good used for RX_Calibrataion channel
        chdata(i).uneq_imp_response=chdata(i).uneq_imp_response*chdata(i).A; % adjust IRx for amplitude
        [chdata(i).uneq_CD_imp_response, ...
            chdata(i).t_DC, ...
            chdata(i).causality_correction_DC_dB, ...
            chdata(i).truncation__DC_dB] = s21_to_impulse_DC(chdata(i).sdc21 ,chdata(i).faxis, param.sample_dt, OP) ;
    end
    % adjust voltage derive here once it's decided what to use
    %------------------------------------------------------------
    % next find Pulse response (SBR) for each channel h^(k)(t)
    if ~OP.DISPLAY_WINDOW && i==1, fprintf('processing COM PDF '); end

    chdata(i).uneq_pulse_response=filter(ones(1, param.samples_per_ui), 1, chdata(i).uneq_imp_response);
    chdata(i).uneq_pulse_DC_response=filter(ones(1, param.samples_per_ui), 1, chdata(i).uneq_CD_imp_response);
    chdata(i).uneq_pulse_CD_response=chdata(i).uneq_pulse_DC_response*chdata(i).A;
    if 1 % not really CD but DC = DC if the channel already has a the Tx added
        % really need to add eq to the DC responce to calc rss.  This is a first pass estimate
        rss=-inf;
        for im=1:param.samples_per_ui
            rss=max(rss, norm( chdata(i).uneq_pulse_CD_response(im:param.samples_per_ui:end)));
        end
        chdata(i).CD_CM_RMS=rss*sqrt(param.sigma_X);
        chdata(i).VCM_HF_struct= get_cm_noise(param.samples_per_ui,chdata(i).uneq_pulse_CD_response,param.levels,param.specBER);
        chdata(i).SCMR=10*log10(max(chdata(1).uneq_pulse_response)^2/chdata(i).VCM_HF_struct.DCn^2);
    end
    if OP.DEBUG && OP.DISPLAY_WINDOW
        if OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
            figure(150+case_number);set(gcf,'Tag','COM');
            screen_size=get(0,'ScreenSize');
            pos = get(gcf, 'OuterPosition');
            set(gcf, 'OuterPosition', ...
                screen_size([3 4 3 4]).*[1 1 0 0] + pos([3 4 3 4]).*[-1 -1 1 1] ...
                - (case_number-1)*[0 20 0 0]);
            %movegui(gcf,'northeast')

            set(gcf, 'Name', sprintf('Case %d PR & PDF - %s', case_number, chdata(i).base));
            subplot(2,1,1);  hold on; % all plots on the same axes
            hp=plot(chdata(i).t, chdata(i).uneq_pulse_response,'Disp', chdata(i).base);
            hold on; % leave on for s-parameter problem finding. RIM 10-02-2023
            hp1=plot(chdata(i).t_DC, chdata(i).uneq_pulse_CD_response,'Disp', [ 'CD ' chdata(i).base ]) ;
        end
        % hide thru PR in order to show xtalk in a reasonable
        % scale. thru is shown in another plot.
        if isequal(chdata(i).type, 'THRU' ) && ~OP.RX_CALIBRATION %|| OP.RX_CALIBRATION % RIM 06-14-2022
            % set(hp, 'visible', 'off');
            % set(get(get(hp,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
        end
        title(sprintf('Unequalized Crosstalk and CD Conversion \n Pulse Responses'))
        ylabel('Volts')
        xlabel('seconds')

        recolor_plots(gca);
    else
        if param.ndfe~=0
            fprintf('%s\tUnequalized pulse peak = %.1f mV\n', chdata(i).base, 1000*max(abs(chdata(i).uneq_pulse_response)));
        end
    end

    fprintf('%s\tCausality correction = %.1f dB', chdata(i).base, chdata(i).causality_correction_dB);
    if OP.ENFORCE_CAUSALITY
        fprintf('\n');
    else
        fprintf(' (not applied)\n');
    end
    fprintf('%s\tTruncation ratio = %.1f dB\n', chdata(i).base, chdata(i).truncation_dB);

end