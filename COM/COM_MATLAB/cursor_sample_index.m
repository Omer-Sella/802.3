function [cursor_i,no_zero_crossing,sbr_peak_i,zxi]=cursor_sample_index(sbr,param,OP,peak_search_range)

%IN:
%sbr = pulse response
%param = COM "param" struct
%OP = COM "OP" struct
%peak_search_range= a limited range to search for the peak (for speed up)
%       it is usually +/- 20 UI
%
%OUT:
%cursor_i = sampling location
%no_zero_crossing = flag that reveals if sampling is not possible.
%       When this function is called in optimize_fom, it signals to quit the current case
%sbr_peak_i = index of the pulse peak.  Note:  tens of thousands of calls to max(sbr) are very
%       time consuming, so saving the peak in one spot is advantageous
%zxi = zero crossing index (returned because RXFFE uses it)

no_zero_crossing=0;
%need to set cursor_i to empty in case no_zero_crossing flag is set
cursor_i=[];

%get peak of pulse and peak index
[max_of_sbr, sbr_peak_tmp]=max(sbr(peak_search_range));
sbr_peak_i=sbr_peak_tmp+peak_search_range(1)-1;


% initial guess at cursor location (t_s)  - based on approximate zero crossing
%limit search space for speed up
search_start=sbr_peak_i-4*param.samples_per_ui;
if search_start<1
    search_start=1;
end
%Find zero crossings
zxi = find(diff(sign(sbr(search_start:sbr_peak_i)-.01*max_of_sbr))>=1)+search_start-1;

%Note:  the original implementation of zxi:
% zxi = find(diff(sign(sbr-.01*max(sbr)))>=1);
% zxi = zxi(zxi<sbr_peak_i);
% zxi = zxi(sbr_peak_i - zxi < 4*param.samples_per_ui);
% The changes to limit search space and remember max(sbr) give 10x speed up
% A test case of 25k runs, reduced from 1.2s to 0.1s


if isempty(zxi)
    %if no zero crossing, the calling program must respond (since sample point will be empty)
    no_zero_crossing=1;
    return;
elseif length(zxi)>1
    %only need the last zero crossing
    zxi=zxi(end);
end
if param.ndfe==0
    max_dfe1=0;
else
    max_dfe1=param.bmax(1);
end
% adjust cursor_i to Solve equation 93A-25 %%
% Muller-Mueller criterion with DFE
mm_range = zxi+(0:2*param.samples_per_ui);
switch OP.CDR
    case 'Mod-MM'
        mm_metric = ...
            abs(sbr(mm_range+param.samples_per_ui)-max_dfe1*sbr(mm_range));
    otherwise % MM
        %MM is generally: first precursor = 0
        %but the actual requirement is for first precursor = first postcursor (after DFE is applied)
        %if first postcursor doesn't exceed max_dfe, then this equates to first precursor = 0
        %in cases where first postcursor exceeds max_dfe, the mismatch is balanced out so that
        %first precursor = first postcursor - max_dfe
        mm_metric = ...
            abs(sbr(mm_range-param.samples_per_ui) - max(sbr(mm_range+param.samples_per_ui)-max_dfe1*sbr(mm_range), 0));
end
[~, mm_cursor_offset] = min(mm_metric);

%cursor_i = the sample location
cursor_i = zxi+mm_cursor_offset-1;