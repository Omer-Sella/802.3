%% floating DFE taps
function [Sout] = interp_Sparam(Sin,fin,fout, ...
    opt_interp_Sparam_mag, opt_interp_Sparam_phase,OP)
% Sout = interp_Sparam(Sin,fin,fout)
%
% Interpolate S-parameters Sin from frequency grid fin to frequency grid
% fout.

if ( fin(end)<fout(end) )
    %    warning('Channel high frequencies extrapolation might be inaccurate!');
end

H_mag = abs(Sin);
H_mag(H_mag<eps)=eps; % handle ill cases...
H_ph = unwrap(angle(Sin));
% For long delay channels, the result can turn anti-causal if frequency step is too coarse. Don't let the
% user ignore that.
if mean(diff(H_ph))>0
    if OP.DEBUG
        warning('Anti-causal response found. Finer frequency step is required for this channel');
    else
        error('Anti-causal response found. Finer frequency step is required for this channel');
    end
end

%opt_interp_Sparam_mag='linear_trend_to_DC';
switch opt_interp_Sparam_mag
    case {'linear_trend_to_DC','linear_trend_to_DC_log_trend_to_inf'}
        if -iscolumn(H_mag), H_mag=H_mag.';end
        if -iscolumn(fin), fin=fin.';end
        fin_x=fin;
        H_mag_x=H_mag(:);
        if fin(1)>0
            p=polyfit(fin(1:10), H_mag(1:10), 1);
            dc_trend_val=polyval(p, 0);
            fin_x=[0, fin_x];
            H_mag_x = [dc_trend_val; H_mag_x];
        end
        if fin(end)<fout(end)
            warn_state=warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
            mid_freq_ind=round(length(fin)/2);
            p=polyfit(fin(mid_freq_ind:end), H_mag(mid_freq_ind:end), 1);
            warning(warn_state);
            hf_trend_val=polyval(p, fout(end));
            if hf_trend_val>H_mag(end)
                hf_trend_val=H_mag(end);
                hf_logtrend_val = H_mag(end);
            elseif hf_trend_val<eps
                hf_trend_val=eps;
                hf_logtrend_val = realmin;
            end
            fin_x=[fin_x, fout(end)];
            H_mag_x = [H_mag_x; hf_trend_val];
        end
        H_mag_i = interp1(fin_x, H_mag_x, fout, 'linear', 'extrap');
        if strcmp(opt_interp_Sparam_mag,'linear_trend_to_DC_log_trend_to_inf')
            H_logmag_i = exp(interp1(fin_x, log([H_mag_x(1:end-1);hf_logtrend_val]), fout, 'linear', 'extrap'));
            indx = find(fout > fin(end),1,'first');
            H_mag_i(indx:end) = H_logmag_i(indx:end);
        end
    case 'trend_to_DC'
        % extrapolate to trend value at DC.
        if -iscolumn(H_mag), H_mag=H_mag.';end
        if -iscolumn(fin), fin=fin.';end
        fin_x=fin;
        H_mag_x=H_mag;
        if fin(1)>0
            warn_state=warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
            p=polyfit(fin(1:10), log10(H_mag(1:10)), 1);
            dc_trend_val=10^polyval(p, 0);
            fin_x=[0, fin_x];
            H_mag_x = [dc_trend_val H_mag_x];
        end
        if fin(end)<fout(end)
            warn_state=warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
            mid_freq_ind=round(length(fin)/2);
            p=polyfit(fin(mid_freq_ind:end), log10(H_mag(mid_freq_ind:end)), 1);
            warning(warn_state);
            hf_trend_val=10^polyval(p, fout(end));
            if hf_trend_val>H_mag(end)
                hf_trend_val=H_mag(end);
            end
            fin_x=[fin_x, fout(end)];
            H_mag_x = [H_mag_x hf_trend_val];
        end
        H_mag_i = 10.^interp1(fin_x,log10(H_mag_x),fout,'linear', 'extrap');
    case 'extrap_to_DC_or_zero'
        % same as extrap_to_DC but detect AC-coupled channels and
        % extrapolate them to 0.
        if fin(1)>0 && 20*log10(H_mag(1))<-20
            % assume AC coupling, with 0 at DC
            H_mag_i = 10.^interp1([0, fin],[-100; log10(H_mag)],fout(fout<=fin(end)),'linear', 'extrap');
        else
            H_mag_i = 10.^interp1(fin, log10(H_mag), fout(fout<=fin(end)),'linear', 'extrap');
        end
        H_mag_i(fout>fin(end)) = H_mag(end);
    case 'extrap_to_DC'
        % first extrapolate down to DC, then use highest available frequency
        % for higher frequencies
        H_mag_i = 10.^interp1(fin,log10(H_mag),fout(fout<=fin(end)),'linear', 'extrap');
        H_mag_i(fout>fin(end)) = H_mag(end);
    case 'old'
        H_mag_i = interp1(fin,H_mag,fout,'linear','extrap');
    otherwise
        error('COM:Extrap:InvalidOption', 'opt_interp_Sparam_mag valid values are "old", "extrap_to_DC"');
end

H_ph_i = interp1(fin,squeeze(H_ph),fout,'linear', 'extrap');

%opt_interp_Sparam_phase='trend_and_shift_to_DC';
%opt_interp_Sparam_phase='interp_cubic_to_dc_linear_to_inf';
switch opt_interp_Sparam_phase
    case 'old'
        H_ph_i = H_ph_i-H_ph_i(1);
    case 'zero_DC'
        H_ph_i(1) = 0;
    case 'interp_to_DC'
        if fin(1) ~= 0
            H_ph_i = interp1([0; fin(:)], [0; H_ph(:)], fout, 'linear', 'extrap');
        end
    case 'extrap_cubic_to_dc_linear_to_inf'
        if fin(1) ~= 0
            % estimate low frequency group delay
            group_delay = -diff(H_ph(:))./diff(fin(:));
            low_freq_gd = group_delay(1:50);
            %  calculate trend, throwing away outliers
            m = median(low_freq_gd); sigma = std(low_freq_gd);
            lf_trend = mean(low_freq_gd(abs(low_freq_gd-m)<sigma));
            % correct outliers in first 10 phase samples
            for k=10:-1:1
                H_ph(k) = H_ph(k+1) + lf_trend*(fin(k+1)-fin(k));
            end
            H_ph_cubic = interp1(fin, H_ph, fout, 'pchip', 'extrap');
            H_ph_linear = interp1(fin, H_ph, fout, 'linear', 'extrap');
            % modification - trend to inf
            if (1)
                high_freq_gd = group_delay(end-50:end);
                %  calculate trend, throwing away outliers
                m = median(high_freq_gd); sigma = std(high_freq_gd);
                hf_trend = -mean(high_freq_gd(abs(high_freq_gd-m)<sigma));
                hf_extrap_range = find(fout>fin(end));
                last_data_sample = hf_extrap_range(1)-1;
                H_ph_linear(hf_extrap_range) = H_ph_linear(last_data_sample) + (fout(hf_extrap_range)-fout(last_data_sample))*hf_trend;
                %                for k=hf_range
                %                    H_ph_linear(k) = H_ph_linear(k-1) + hf_trend*(fout(k)-fout(k-1));
                %                end
            end
            
            [UNUSED_OUTPUT, indx] = min(abs(H_ph_cubic-H_ph_linear)); %#ok<ASGLU>
            H_ph_i=H_ph_cubic;
            H_ph_i(indx:end) = H_ph_linear(indx:end);
            H_ph_i = H_ph_linear; % John Ewen 12/13/2019
        end
    case 'interp_and_shift_to_DC'
        if fin(1) ~= 0
            dc_phase_trend = H_ph(1)-(H_ph(2)-H_ph(1))/(fin(2)-fin(1))*fin(1);
            H_ph_i = interp1([0; fin(:)], [0; H_ph(:)-dc_phase_trend], fout, 'linear', 'extrap');
        end
    case 'trend_and_shift_to_DC'
        % estimate low frequency group delay
        group_delay = -diff(H_ph(:))./diff(fin(:));
        low_freq_gd = group_delay(1:50);
        %  calculate trend, throwing away outliers
        m = median(low_freq_gd); sigma = std(low_freq_gd);
        lf_trend = mean(low_freq_gd(abs(low_freq_gd-m)<sigma));
        fin_x=fin;
        H_ph_x=H_ph(:);
        if fin(1) ~= 0
            % correct outliers in first 10 phase samples
            for k=10:-1:1
                H_ph(k) = H_ph(k+1) + lf_trend*(fin(k+1)-fin(k));
            end
            
            % shift all phase data so that DC extrapolation to 0 follows trend
            dc_phase_trend = H_ph(1)+lf_trend*(fin(1)-0);
            fin_x=[0, fin_x];
            H_ph_x=[0; H_ph(:)-dc_phase_trend];
        end
        % Modification: extrapolate using trend. (interp1 with "extrap" extrapolates using just
        % the last two samples, so noise can create an inverted slope and
        % non-causal response).
        if fout(end)>fin(end)
            group_delay = -diff(H_ph_x(:))./diff(fin_x(:));
            % p=polyfit(fin_x', H_ph_x, 1);
            hf_phase_trend = H_ph_x(end)-median(group_delay)*(max(fout)-max(fin_x));
            % hf_phase_trend=polyval(p,max(fout));
            fin_x=[fin_x, fout(end)];
            H_ph_x=[H_ph_x; hf_phase_trend];
        end
        H_ph_i = interp1(fin_x, H_ph_x, fout, 'linear', 'extrap');
        
    otherwise
        error('COM:Extrap:InvalidOption', ...
            'debug_interp_Sparam valid values are "old", "zero_DC", "interp_to_DC", "interp_and_shift_to_DC", "trend_and_shift_to_DC", "interp_cubic_to_dc_linear_to_inf"');
end
H_i = H_mag_i.*exp(1j*H_ph_i);
Sout=H_i;