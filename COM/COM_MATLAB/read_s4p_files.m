function [chdata,SDDch,SDDp2p] = read_s4p_files(param, OP, chdata)
%% extract s-parameter and convert to differential mode
% extract s-parameter data from files and apply tx and rx filters as well as package filters
num_files=length(chdata);
if ~OP.DISPLAY_WINDOW, fprintf('reading file '); end
for i=1:num_files
    if OP.DISPLAY_WINDOW; hwaitbar=waitbar(0);end
    progress = i/num_files;
    if OP.DISPLAY_WINDOW
        [~,a]=fileparts(chdata(i).filename);
        waitbar(progress, hwaitbar, ['Processing ' a]); figure(hwaitbar); drawnow;
    else
        fprintf('%i ',i);
    end
    
    % Skip reading file if it was already read (multiple test cases)
    if (~isfield(chdata(i), 'faxis')) || isempty(chdata(i).faxis)
        switch lower(chdata(i).ext)
            case '.s2p' % for differential return loss
                [Sch,SDDch] = read_p2_s2params(chdata(i).filename,  0, 0, param.snpPortsOrder, OP);
                chdata(i).fmaxi = length(Sch.freq);
                chdata(i).faxis = Sch.freq;
                chdata(i).sdd11_raw = transpose(SDDch(1:chdata(i).fmaxi,1,1));
                SDDp2p(i)=NaN;
                chdata(i).sdd11_orig=chdata(i).sdd11_raw;
                chdata(i).sdd11=chdata(i).sdd11_raw;
            case '.s4p'
                if length(param.snpPortsOrder) ~= 4
                    error( 'warning:sNpFilePortMismatch', ...
                        '\n\t The number of ports defined (%G) does not match the sNp file type (%s)', ...
                        length(param.snpPortsOrder), ...
                        chdata(i).ext ...
                        );
                end
                % read function returns differnetial mode parameters
                if param.package_testcase_i==1 % added to speed up cases e.g. don't read file in twice
                    [Sch, SDDch, SDCch] = read_p4_s4params(chdata(i).filename,  0, 0, param.snpPortsOrder, OP,param);
%                     param.holdsdata(i).Sch= Sch;
%                     param.holdsdata(i).SDDch=  SDDch;
%                     param.holdsdata(i).SDCch= SDCch;
                else
                    error('If this line is reached, there is a logic error');
%                     Sch=param.holdsdata(i).Sch;
%                     SDDch=param.holdsdata(i).SDDch;
%                     SDCch=param.holdsdata(i).SDCch;
                end
                chdata(i).fmaxi = length(Sch.freq);
                

                if Sch.freq(chdata(i).fmaxi) < param.fb
                    warning('COM:read_s4p:MaxFreqTooLow', ...
                        'In %s: the maximum frequency provided, %g, is less than the signaling rate: %g', ...
                        chdata(i).filename, Sch.freq(end), param.fb);
                end
                if Sch.freq(1) > param.max_start_freq
                    warning('COM:read_s4p:StartFreqTooHigh', ...
                        'In %s: minimum frequency, %.2g GHz, is larger than the recommended %.2g GHz', ...
                        chdata(i).filename, Sch.freq(1)/1e9, param.max_start_freq/1e9);
                end
                freqstep=diff(Sch.freq);
                % ignore frequency differences up to 1 Hz - possible numerical artifacts
                if max(freqstep)-min(freqstep) > 1
                    warning('COM:read_s4p:NonUniformFreqSpacing', 'In %s: non-uniform frequency steps: min=%.3g GHz, max=%.3g GHz', ...
                        chdata(i).filename, min(freqstep)/1e9, max(freqstep)/1e9);
                end
                if max(freqstep) - param.max_freq_step > 1
                    warning('COM:read_s4p:FreqStepTooHigh', 'In %s: frequency step, %.2g GHz, is larger than the recommended %.2g GHz', ...
                        chdata(i).filename, max(freqstep)/1e9, param.max_freq_step/1e9);
                end
                
                chdata(i).faxis = Sch.freq;
                chdata(i).sdd12_raw = transpose(SDDch(1:chdata(i).fmaxi,1,2));
                chdata(i).sdd21_raw = transpose(SDDch(1:chdata(i).fmaxi,2,1));
                chdata(i).sdd22_raw = transpose(SDDch(1:chdata(i).fmaxi,2,2));
                chdata(i).sdd11_raw = transpose(SDDch(1:chdata(i).fmaxi,1,1));
                % mode conversion
                chdata(i).sdc12_raw = transpose(SDCch(1:chdata(i).fmaxi,1,2));
                chdata(i).sdc21_raw = transpose(SDCch(1:chdata(i).fmaxi,2,1));
                chdata(i).sdc22_raw = transpose(SDCch(1:chdata(i).fmaxi,2,2));
                chdata(i).sdc11_raw = transpose(SDCch(1:chdata(i).fmaxi,1,1));
                %save original and add board (if required)
                chdata(i).sdd11_orig=chdata(i).sdd11_raw;
                chdata(i).sdd22_orig=chdata(i).sdd22_raw;
                chdata(i).sdd12_orig=chdata(i).sdd12_raw;
                chdata(i).sdd21_orig=chdata(i).sdd21_raw;
                if OP.include_pcb
                    % add boards to sdd
                    [chdata(i).sdd11_raw, chdata(i).sdd12_raw, chdata(i).sdd21_raw, chdata(i).sdd22_raw] = add_brd(chdata(i), param, OP);
                    
                end
                %save final return loss (after the boards were included)
                chdata(i).sdd11=chdata(i).sdd11_raw;
                chdata(i).sdd22=chdata(i).sdd22_raw;
            otherwise
                error('Extension "%s" in file "%s" is not supported',chdata(i).ext,chdata(i).filename);
        end
        
        %Crosstalk frequency axis must be the same as Thru
        if i>1
            %error on length difference
            if length(chdata(i).faxis)~=length(chdata(1).faxis)
                error('Crosstalk file "%s" has different number of frequency points',chdata(i).filename);
            end
            %error if any value > 1Hz (don't want to check for exact
            %equality in case of floating point error)
            Fdiff=abs(chdata(i).faxis-chdata(1).faxis);
            if max(Fdiff)>1
                error('Crosstalk file "%s" has a different frequency axis',chdata(i).filename);
            end
        end
    else
        SDDch(:,1,2)=chdata(i).sdd12_raw;
        SDDch(:,2,1)=chdata(i).sdd21_raw;
        SDDch(:,1,1)=chdata(i).sdd11_raw;
        SDDch(:,2,2)=chdata(i).sdd22_raw;
    end
    chdata(i).sigma_ACCM_at_tp0=0;
    if ~param.FLAG.S2P
        if  OP.INC_PACKAGE ~= 0 || (OP.RX_CALIBRATION == 1  && i==1)
            if (OP.RX_CALIBRATION == 1  && i==2)
                chdata(i).sdd21=chdata(i).sdd21_raw;
            else
                %updated package construction with single function for both DD and DC
                [chdata(i).sdd21p,SDDp2p(i)]= s21_pkg(chdata(i), param, OP, i);
                [chdata(i).sdd21p_nodie]= s21_pkg(chdata(i), param, OP, i, 'dd', 0);
                chdata(i).sdd21=chdata(i).sdd21p;
                if 1 % for AC CM noise inclusion
                    [chdata(i).sdc21p,SDCp2p(i),chdata(i).sigma_ACCM_at_tp0]= s21_pkg(chdata(i), param, OP, i,'dc');
                    chdata(i).sdc21=chdata(i).sdc21p;
                end
            end
        else
            chdata(i).sdd21=chdata(i).sdd21_raw;
        end
        chdata(i).sdd21f=chdata(i).sdd21_orig; % used for FD analysis i.e. not filtered (RIM 9/24/2021 without boards or packages)
    end
end
if ~OP.DISPLAY_WINDOW, fprintf('\n'); end
