function [p_burst,p_error_propagation]=Burst_Probability_Calc(COM_SNR_Struct,DFE_taps,param,OP)

% an error burst of length N will cause each of the first N taps tap to mis-correct and create a PAM (2 or
% 4) noise term - depending on the N'th previous symbol, with double the tap voltage. From this we calculate
% the probability of staying in error state, i.e. burst of length N+1.

A_s=COM_SNR_Struct.A_s;
% initialize loop with uncorrelated noise and BER
error_propagation_noise_pdf{1}=COM_SNR_Struct.combined_interference_and_noise_pdf;  % PDF for burst of length 1 is the uncorrelated PDF

% Assume an error will occur if the noise excceds the available signal
% reduced by some dB. reduction is by COM threshold minus Error
% propagation margin (a positive EP margin reduces uncorrelated error probability
% below target BER).
error_threshold =  A_s./10^((param.pass_threshold-OP.COM_EP_margin)/20);
% Find the probability of this event by integration of the PDF. Use 1e-20 as a floor probabilty if noise PDF isn't wide enough.
x_error_propagation = find(error_propagation_noise_pdf{1}.x >= error_threshold, 1, 'first');
if isempty(x_error_propagation)
    p_error_propagation(1) = 1e-20;
else
    p_error_propagation(1) = sum(error_propagation_noise_pdf{1}.y(x_error_propagation:end));  % uncorrelated BER
end

sorted_abs_dfe_taps = sort(abs(DFE_taps), 'descend');
for k=2:min(param.ndfe, OP.nburst)
    % (arrays kept to allow tracking during development, though not really needed)
    if OP.use_simple_EP_model
        post_error_dfe_noise_pdf{k} = get_pdf_from_sampled_signal( 2*A_s*max(sorted_abs_dfe_taps), param.levels, param.delta_y ); %#ok<AGROW>
        error_propagation_noise_pdf{k} = conv_fct(error_propagation_noise_pdf{1}, post_error_dfe_noise_pdf{k}); %#ok<AGROW>
    else
        post_error_dfe_noise_pdf{k} = get_pdf_from_sampled_signal( 2*A_s*sorted_abs_dfe_taps(k-1), param.levels, param.delta_y ); %#ok<AGROW>
        error_propagation_noise_pdf{k} = conv_fct(error_propagation_noise_pdf{k-1}, post_error_dfe_noise_pdf{k}); %#ok<AGROW>
    end
    
    % Assume an error will propagate if this noise exceeds the threshold defined above
    x_error_propagation = find(error_propagation_noise_pdf{k}.x >= error_threshold, 1, 'first');
    if isempty(x_error_propagation)
        p_error_propagation(k) = 1e-20; %#ok<AGROW>
    else
        p_error_propagation(k) = sum(error_propagation_noise_pdf{k}.y(x_error_propagation:end)); %#ok<AGROW>
    end
end

% Assume an uncorrelated error will occur if the original noise exceeds
% the available signal reduced by pass_threhsold dB. Find the probability
% of this event by partial sum of the PDF.
% p_uncorrelated_error_i = find(combined_interference_and_noise_pdf.x >= A_s./10^(param.pass_threshold/20), 1, 'first');
% p_uncorrelated_error = sum(combined_interference_and_noise_pdf.y(p_uncorrelated_error_i:end));

% probability of bursts of different lengths
p_burst = cumprod(p_error_propagation);