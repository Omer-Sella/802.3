function [msg] = end_display_control(msg,param,OP,output_args,COM,min_ERL,ERL,VEO_mV,VEC_dB,threshold_DER,DISPLAY_WINDOW)
[ncases, mele]=size(param.z_p_tx_cases);
if mele ==2
    param.flex=2;
elseif mele==4
    param.flex=4;
elseif mele==1
    param.flex=1;
else
    error(springf('config file syntax error'))
end


if DISPLAY_WINDOW && ~OP.RX_CALIBRATION
    % display bathtub curves in one axis per test case.
    %             h=findall(0, 'Name', 'COM results');
    if ~exist('h','var')
        msgtext = cell(1, length(OP.pkg_len_select));
        msgcolor = 'g';
    else
        msgtext=get(findobj(h, 'type', 'text'), 'string');
        msgcolor = get(h, 'color');
        close(h); % will be recreated
    end
    msgctr=size(msgtext,1)+1;
    if ~OP.ERL_ONLY
        switch OP.PHY
            case 'C2M'
                if VEO_mV >= param.Min_VEO && VEO_mV <= param.Max_VEO
                    msg=sprintf('%s: EH = %.3f mV (pass)\n', ...
                        msg, VEO_mV);
                else
                    msg=sprintf('%s: EH = %.3f mV (FAIL)\n', ...
                        msg, VEO_mV);
                    msgcolor = 'r';
                end
                
                if VEC_dB <= param.VEC_pass_threshold
                    msg=sprintf('%s: VEC = %.3f dB (pass)\n', ...
                        (msg), VEC_dB);
                else
                    msg=sprintf('%s: VEC = %.3f dB (FAIL)\n', ...
                        (msg), VEC_dB);
                    msgcolor = 'r';
                end
            case 'C2C'
                if COM >= param.pass_threshold
                    %                         msgtext{package_testcase_i}=sprintf('%s: COM = %.3f dB (pass)\n', ...
                    %                             msg, COM);
                    msg=sprintf('%s: COM = %.3f dB (pass)\n', ...
                        msg, COM);
                else
                    %                         msgtext{package_testcase_i}=sprintf('%s: COM = %.3f dB (FAIL)\n', ...
                    %                             msg, COM);
                    msg=sprintf('%s: COM = %.3f dB (FAIL)\n', ...
                        msg, COM);
                    msgcolor = 'r';
                end
                % begin yasuo patch 3/18/2019
                msg=sprintf('%s: DER = %.3e at COM threshold \n', msg, threshold_DER);
                % end yasuo patch
            case 'C2Mcom'
                if VEO_mV >= param.Min_VEO && VEO_mV <= param.Max_VEO
                    msg=sprintf('%s: EH = %.3f mV (pass)\n', ...
                        msg, VEO_mV);
                else
                    msg=sprintf('%s: EH = %.3f mV (FAIL)\n', ...
                        msg, VEO_mV);
                    msgcolor = 'r';
                end
                
                if VEC_dB <= param.VEC_pass_threshold
                    msg=sprintf('%s: VEC = %.3f dB (pass)\n', ...
                        (msg), VEC_dB);
                else
                    msg=sprintf('%s: VEC = %.3f dB (FAIL)\n', ...
                        (msg), VEC_dB);
                    msgcolor = 'r';
                end
                if COM >= param.pass_threshold
                    %                         msgtext{package_testcase_i}=sprintf('%s: COM = %.3f dB (pass)\n', ...
                    %                             msg, COM);
                    msg=sprintf('%s: COM = %.3f dB (pass)\n', ...
                        msg, COM);
                else
                    %                         msgtext{package_testcase_i}=sprintf('%s: COM = %.3f dB (FAIL)\n', ...
                    %                             msg, COM);
                    msg=sprintf('%s: COM = %.3f dB (FAIL)\n', ...
                        msg, COM);
                    msgcolor = 'r';
                end
                % begin yasuo patch 3/18/2019
                msg=sprintf('%s: DER = %.3e at COM threshold \n', msg, threshold_DER);
                % end yasuo patch
        end
    end
    if OP.ERL
        if ~isempty(ERL)
            if min_ERL >= param.ERL_pass_threshold
                %                             msgtext{package_testcase_i+1}=[sprintf(' PASS ... ERL = %.3f dB\n', min_ERL)];
                msg=[sprintf('%s: PASS ... ERL = %.3f dB (%.3f dB,%.3f dB) \n',msg, min_ERL, ERL)];
            else
                %                             msgtext{package_testcase_i+1}=[ sprintf(' FAIL ... ERL = %.3f dB\n', min_ERL) ];
                msg=[ sprintf('%s: FAIL ... ERL = %.3f dB (%.3f dB,%.3f dB) \n',msg, min_ERL, ERL) ];
                msgcolor = 'r';
            end
        end
    end
    h=msgbox(msg, ['COM r' output_args.code_revision ' results']);
    set(h, 'color', msgcolor, 'tag', 'COM');
    movegui(h, 'center');
else % no windows
    %     display(['max noise at BER = ' num2str(peak_interference_at_BER)])
    %     display(['signal after eq = ' num2str(A_s/(param.levels-1))])
    if ~OP.ERL_ONLY
        switch OP.PHY
            case 'C2C'
                if COM >= param.pass_threshold
                    fprintf('%s <strong> PASS ... COM = %.3f dB</strong>\n', msg, COM);
                else
                    fprintf(2,'%s <strong> FAIL ... COM = %.3f dB</strong>\n',  msg, COM);
                end
                % begin yasuo patch 3/18/2019
                fprintf('%s DER = %.3e at COM threshold \n', msg, threshold_DER);
                % end yasuo patch
            case 'C2Mcom'
                if VEC_dB <= param.VEC_pass_threshold
                    fprintf('%s <strong> PASS ... VEC = %.3f dB</strong>\n', msg, VEC_dB);
                else
                    fprintf(2,'%s <strong> FAIL ... VEC = %.3f dB</strong>\n',  msg, VEC_dB);
                end
                
                if VEO_mV >= param.Min_VEO && VEO_mV <= param.Max_VEO
                    fprintf('%s <strong> PASS ... EH = %.3f mV</strong>\n', msg, VEO_mV);
                else
                    fprintf(2,'%s <strong> FAIL ... EH = %.3f mV</strong>\n',  msg, VEO_mV);
                end
                if COM >= param.pass_threshold
                    fprintf('%s <strong> PASS ... COM = %.3f dB</strong>\n', msg, COM);
                else
                    fprintf(2,'%s <strong> FAIL ... COM = %.3f dB</strong>\n',  msg, COM);
                end
                % begin yasuo patch 3/18/2019
                fprintf('%s DER = %.3e at COM threshold \n', msg, threshold_DER);
                % end yasuo patch
            case 'C2M'
                if VEC_dB <= param.VEC_pass_threshold
                    fprintf('%s <strong> PASS ... VEC = %.3f dB</strong>\n', msg, VEC_dB);
                else
                    fprintf(2,'%s <strong> FAIL ... VEC = %.3f dB</strong>\n',  msg, VEC_dB);
                end
                
                if VEO_mV >= param.Min_VEO && VEO_mV <= param.Max_VEO
                    fprintf('%s <strong> PASS ... EH = %.3f mV</strong>\n', msg, VEO_mV);
                else
                    fprintf(2,'%s <strong> FAIL ... EH = %.3f mV</strong>\n',  msg, VEO_mV);
                end
        end
    end
    if OP.ERL
        if ~isempty(ERL)
            if min_ERL >= param.ERL_pass_threshold
                %                             msgtext{package_testcase_i+1}=[sprintf(' PASS ... ERL = %.3f dB\n', min_ERL)];
                fprintf('%s: <strong> PASS ... ERL = %.3f dB (%.3f dB, %.3f dB)</strong>\n',msg, min_ERL, ERL );
            else
                %                             msgtext{package_testcase_i+1}=[ sprintf(' FAIL ... ERL = %.3f dB\n', min_ERL) ];
                fprintf(2,'%s:  <strong> FAIL ... ERL = %.3f dB (%.3f dB, %.3f dB)</strong>\n',msg, min_ERL, ERL) ;
            end
        end
    end
end
