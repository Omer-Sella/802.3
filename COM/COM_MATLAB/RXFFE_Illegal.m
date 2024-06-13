function is_illegal=RXFFE_Illegal(C,param,last_index)

%check if RXFFE taps are illegal
%C = RXFFE taps
%param = COM param struct
%last_index is used when computing illegality prior to Backoff.  It will be set so taps
% in the Backoff region are not considered in the legality check.

%If last index is omitted, set it to length(C)
if nargin<3
    last_index=length(C);
end

is_illegal=0;

%Check cursor tap
Ccur_i=param.RxFFE_cmx+1;
if C(Ccur_i) < param.ffe_main_cursor_min
    is_illegal=1;
    return;
end

%Check postcursors
if param.ffe_post_tap_len ~=0
    if abs(C(Ccur_i +1)) > param.ffe_post_tap1_max
        is_illegal=1;
        return;
    end
    if (param.ffe_post_tap_len > 1)
        if sum(abs(C((Ccur_i +2):last_index))  > param.ffe_tapn_max)
            is_illegal=1;
            return;
        end
    end
end

%Check precursors
if param.ffe_pre_tap_len ~=0
    if abs(C(Ccur_i -1)) > param.ffe_pre_tap1_max
        is_illegal=1;
        return;
    end
    if (param.ffe_pre_tap_len > 1)
        %                                             if sum(abs(C((Ccur_i +2):end)) > param.ffe_tapn_max) , continue; end
        if sum(abs(C(1:(Ccur_i - 2))) > param.ffe_tapn_max)
            is_illegal=1;
            return;
        end % 11.22.2018 Yasou Hadaka
    end
end