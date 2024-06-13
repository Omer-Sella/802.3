function [delay_sec, delay_idx]= calculate_delay_CausalityEnforcement(freq, sdd21, param, OP)
% History:
% 1. 14th October, 2021 (Intial release)

% Definition:
% This function captures the channel delay through the time domain using causality enforcement. 
% Following are the steps being followed.
% Step 1. Cascade negative frequencies
% Step 2. Extract magnitude
% Step 3. IFFT of the magnitude
% Step 4. Multiply by the sign(t)
% Step 5. Calculate the phase of the 1j*causal_phase
% Step 6. casual_function= |original|*exp(-1j*causal_phase)
% Step 7. f-domain to t-domain pulse response
% Step 8. Calculate the delay

% Author:
% Hansel Dsilva (dsilvahansel@gmail.com or hanseldsilva@achronix.com)

% Reference:
% 1] "IEEE Standard for Electrical Characterization of Printed Circuit Board and Related Interconnects at Frequencies up to 50 GHz," in IEEE Std 370-2020 , vol., no., pp.1-147, 8 Jan. 2021, doi: 10.1109/IEEESTD.2021.9316329.

% Input:
% freq          %frequency in hertz (odd number points)
% sdd21          %insertion loss in complex (odd number points)
% param          %COM native structure passed
% OP          %COM native structure passed

% Output:
% delay_sec          %channel delay in seconds
% delay_idx          %channel delay index

if iscolumn(sdd21)
    sdd21= sdd21.';
end
if iscolumn(freq)
    freq= freq.';
end

%---start. Step 1. Cascade negative frequencies
% sdd21_conj= zeros(1, length(freq)+ length(freq) - 1);
sdd21_conj= [sdd21, conj(sdd21(end:-1:2))];%For some reason this only works 
% sdd21_conj = [real(sdd21(1)), sdd21(2:end-1), real(sdd21(end)), flipud(conj(sdd21(2:end-1)))];
%---end. Step 1. Cascade negative frequencies

%---start. Step 2. Extract magnitude
sdd21_mag_conj = real(log(abs(sdd21_conj)));
%---end. Step 2. Extract magnitude

%---start. Step 3. IFFT of the magnitude
sdd21_mag_time = ifft(sdd21_mag_conj);
%---end. Step 3. IFFT of the magnitude

%---start. Step 4. Multiply by the sign(t)
sdd21_mag_time(1:length(freq))= +1j*sdd21_mag_time(1:length(freq));
sdd21_mag_time(length(freq)+1:end)= -1j*sdd21_mag_time(length(freq)+1:end);
%---end. Step 4. Multiply by the sign(t)

%---start. Step 5. Calculate the phase of the original_signal and the 1j*causal_phase
sdd21_phase_causality_enforced = real(fft(sdd21_mag_time));
%---end. Step 5. Calculate the phase of the original_signal and the 1j*causal_phase

%---start. Step 6. casual_function= |original|*exp(-1j*causal_phase)
sdd21_causality_enforced= abs(sdd21_conj).*exp(-1j*sdd21_phase_causality_enforced);
sdd21_causality_enforced= sdd21_causality_enforced(1:length(freq));
%---end. Step 6. casual_function= |original|*exp(-1j*causal_phase)

%---start. Step 7. f-domain to t-domain pulse response
%------Note. Do not use s21_to_impulse() for we do not want to truncate
%--------- Extrapolation has been already done by the COM tool
freq_array= freq;
time_step= param.sample_dt;
fmax=1/time_step/2;
freq_step=(freq_array(3)-freq_array(2))/1;
fout=0:1/round(fmax/freq_step)*fmax:fmax;

ILin=sdd21;
IL=interp_Sparam(ILin,freq_array,fout, ...
    OP.interp_sparam_mag, OP.interp_sparam_phase,OP);
IL_nan = find(isnan(IL));
for in=IL_nan
    IL(in)=IL(in-1);
end
IL = IL(:);
% add padding for time steps
IL_symmetric = [real(IL(1)); IL(2:end-1); real(IL(end)); flipud(conj(IL(2:end-1)))];
sdd21_PR = filter(ones(1, param.samples_per_ui), 1, real(ifft(IL_symmetric)));%real(ifft(IL_symmetric));
clear IL IL_nan IL_symmetric

ILin=sdd21_causality_enforced;
IL=interp_Sparam(ILin,freq_array,fout, ...
    OP.interp_sparam_mag, OP.interp_sparam_phase,OP);
IL_nan = find(isnan(IL));
for in=IL_nan
    IL(in)=IL(in-1);
end
IL = IL(:);
% add padding for time steps
% IL_symmetric = [IL(1:end-1); 0; flipud(conj(IL(2:end-1)))];
% IL_symmetric = [IL(1:end); flipud(conj(IL(2:end-1)))];
IL_symmetric = [real(IL(1)); IL(2:end-1); real(IL(end)); flipud(conj(IL(2:end-1)))];
sdd21_causality_enforced_PR = filter(ones(1, param.samples_per_ui), 1, real(ifft(IL_symmetric)));%real(ifft(IL_symmetric));
clear IL IL_nan IL_symmetric

clear time_step fmax freq_step freq_array

freq_step=(fout(3)-fout(2))/1;
L= length(sdd21_PR);
t_base = (0:L-1)/(freq_step*L);
clear fout freq_step L
%---end. Step 7. f-domain to t-domain pulse response

%---start. Step 8. Calculate the delay
%------start. Remove the last 5% of the waveform for noise due to IFFT
sdd21_causality_enforced_PR_reduced= sdd21_PR(1:floor(end*95/100));
sdd21_PR_reduced= sdd21_causality_enforced_PR(1:floor(end*95/100));
%------end. Remove the last 5% of the waveform for noise due to IFFT

%---start. calculate the difference in index between the peaks
[~, peak_x_idx] = max(sdd21_causality_enforced_PR_reduced);
[~, peak_y_idx] = max(sdd21_PR_reduced);
peak_idx_difference = peak_x_idx - peak_y_idx;
%---end. calculate the difference in index between the peaks

if peak_idx_difference~=0
    search_bounds = min(peak_x_idx, peak_y_idx);%<----one may fix it to 1e3; using minimum of the peaks in assuming the other signal has zero delay
    error_value = length(sdd21_causality_enforced_PR_reduced);
    error_idx = 0;
    %     i= 1;
    for shift_value= peak_idx_difference-search_bounds:peak_idx_difference+search_bounds
        sdd21_PR_reduced_shifted = circshift(sdd21_PR_reduced,shift_value);
        current_error = norm(sdd21_PR_reduced_shifted-sdd21_PR_reduced);
        if (error_value > current_error)
            error_idx = shift_value;
            error_value = current_error;
        end
        %         error_idx_H(i)= error_idx;
        %         i= i+ 1;
    end
    %plot(error_idx_H);
    
    if error_idx==0
        error('An odd case has occured when calcualting the channel delay. Please contact the tool developer.');
    end
    
    delay_idx = error_idx;
else
    delay_idx= 1;
end

delay_sec= t_base(abs(delay_idx));