% function [ bmax floating_tap_locations] = floating_taps( hisi,   N_b,     N_bf,  N_bg,  N_bmax,       bmaxg. COO))
function [ Vfiltered, Cmod, idx ]= force( V ,param, OP , ix, C, return_V, chdata, txffe, Noise_XC)
% Vfilter is vector forced filtered sbr
% Cmod is the ffe tap co-efficient vector
% if C is passed, just process V with C else compute C
% cmx=param.rx_cmx; number of pre cursor taps
% cpx=param.rx_cps; number of post cursor taps
% V=sbr; pass pulse response
% ix the sample point in the passed pulse response
% the sample point is recomputed by optimize_fom
% idx - return floating tap location (RIM 9-19-2023)
% OP not used for now
%return_V is a flag with default value = 1.  If 0, Vfiltered is not returned
%   this allows significant speed up in optimize_fom since FFE is time consuming
%   and many combinatiFons of "Cmod" result in illegal combinations that do not need
%   Vfiltered to be calculated
% test with load('SBR_FIR_resp.mydata','-mat')
idx=[];
if nargin<4
    ix=find(V==max(V),1,'first');
end
if nargin<5
    C=[];
end
if nargin<6
    return_V=1;
end
cmx=param.RxFFE_cmx;
cpx=param.RxFFE_cpx;
% do this early on so we can reuse the old code
if param.N_bg ~=0 % must be floating taps
    cpx=param.N_bmax; % N_f in spreadsheet
end
num_taps=cmx+cpx+1;
cstep=param.RxFFE_stepz;
ndfe=param.ndfe;
spui=param.samples_per_ui;
param.current_ffegain=0;
if return_V && ~isempty(C)
    %    RIM 2-3-23 when we just want to EQ not find EQ
    Vfiltered=FFE( C , param.RxFFE_cmx,spui, V );
    Cmod=C;
    return
end
% Aling V to ix ( the sample point) and then create the sampled vector vsampled_raw
if ix < length(V)
    if isrow(V)
        if mod(ix,spui) == 0
            vsampled_raw = [V(spui+mod(ix,spui):spui:(mod(ix,spui)+spui*(floor(ix/spui)-1)))'; V(ix:spui:end)'];
        else
            vsampled_raw = [V(mod(ix,spui):spui:(mod(ix,spui)+spui*(floor(ix/spui)-1)))'; V(ix:spui:end)'];
        end

    else
        if mod(ix,spui) == 0
            vsampled_raw = [V(spui+mod(ix,spui):spui:(mod(ix,spui)+spui*(floor(ix/spui)-1))); V(ix:spui:end)];
        else
            vsampled_raw = [V(mod(ix,spui):spui:(mod(ix,spui)+spui*(floor(ix/spui)-1))); V(ix:spui:end)];
        end
    end
else
    if isrow(V)
        if mod(ix,spui) == 0%Yasou Hidaka 11/16/2018
            vsampled_raw = V(spui+mod(ix,spui):spui:end)';%Yasou Hidaka 11/16/2018
        else
            vsampled_raw = V(mod(ix,spui):spui:end)';%Yasou Hidaka 11/16/2018
        end
    else
        if mod(ix,spui) == 0%Yasou Hidaka 11/16/2018
            vsampled_raw = V(spui+mod(ix,spui):spui:end) ;
        else
            vsampled_raw = V(mod(ix,spui):spui:end) ;%Yasou Hidaka 11/16/2018
        end
    end
end
% zero pad vsampled to account for PR with short delay. RIM 10-02-2023
vsampled=[zeros(1,num_taps) vsampled_raw' zeros(1,cpx)];% pad for pre and post cursor prior to shifting

%% find the index, ivs, for the sample point but in the UI resample vector, vsampled
% ivs=find(vsampled==V(ix),1,'first');% ivs is the sample point for V
% Upen Kareti suggested fix for indexing 11/04/18
if ix < length(V)
    ivs=find(vsampled==V(ix),1,'first');% ivs is the sample point for V
else
    ivs=find(vsampled == max(vsampled),1,'first');
end


%% create VV matrix of shifted UI spaced sample of the pulse response
% only consider the VV matrix that correstonds to the FFE taps
VV=zeros(num_taps,num_taps);
for i=1:num_taps
    start_idx=ivs+i-1;
    end_idx=start_idx-num_taps+1;
    VV(:,i)=vsampled(start_idx:-1:end_idx);
end
% may want to test VV*VV' for rcond here not sure what value to use but 1e-5 is always bad
%% Apply RXFFE
if isempty(C)
    switch upper(OP.FFE_OPT_METHOD)
        case 'WIENER-HOPF'
            C= WIENER_HOPF_MMSE(vsampled ,param, OP , chdata, txffe, Noise_XC,ivs) ;
            Cmod=C(1:num_taps);
        otherwise
            % cmx+1 is the cursor or sample point
            %VV=VV(:,ivs-cmx:ivs+cpx); % only consider the VV matrix that correstonds to the FFE taps
            FV=zeros(1,num_taps); % zero the forceing vector, FV first
            FV(cmx+1)=vsampled(ivs)*10^(param.current_ffegain/20); % force the voltage at sample point
            if param.ndfe~=0 && (cpx > 0) % Yasuo Hidaka suggest fix for no postC 11/11/18
                %          FV(cmx+2)=param.bmax(1)*FV(cmx+1);   % force the post cursor to bmax if dfe exists
                FV(cmx+2)=min(param.bmax(1)*FV(cmx+1),abs(vsampled(ivs+1)))*sign(vsampled(ivs+1));
            end
            %C=((VV'*VV)^-1*VV')'*FV'; % sikve for FFE taps, C
            if diff(size(VV))==0
                %For square matrix, can solve C using simple inv(VV')*FV'
                C=VV'\FV';
            else
                %otherwise use the general solution with psuedo inverse
                %note:  this is the same as doing pinv(VV') but pinv is far slower
                % C=(inv(VV'*VV)*VV')'*FV'; % sikve for FFE taps, C
                C=(inv(VV*VV')*VV)*FV'; % sikve for FFE taps, C
            end

            Cmod=C(1:num_taps);
    end


    % added for 4.2 find floating taps with either ISI or taps
    switch lower(OP.RXFFE_FLOAT_CTL)
        case 'taps'
            [idx]=findbankloc(Cmod,param.N_tail_start,param.N_bmax,param.N_bf,Cmod(cmx+1),param.bmaxg,param.N_bg);
        otherwise
            [idx]=findbankloc(VV(cmx+1,:),param.N_tail_start,param.N_bmax,param.N_bf,Cmod(cmx+1),param.bmaxg,param.N_bg);
    end
    switch lower(OP.RXFFE_TAP_CONSTRAINT)
        case 'unity cursor'
            Cmod=Cmod/Cmod(cmx+1);
        otherwise
            Cmod=C;
    end
    if cstep ~= 0
        Cmod=floor(abs(Cmod/cstep)).*sign(Cmod)*cstep;% r250 quantize with floor ad sign(C)
    end

    if ~isempty(idx)
        idx=sort(idx);
        C1=Cmod;
        % C1(param.N_tail_start:end)=0;
        C1(param.RxFFE_cmx+param.RxFFE_cpx+2:end)=0; % zero out flosting taps
        C1(cmx+1+idx)=Cmod(cmx+1+idx);
        Cmod=C1;
    else
        % Cmod=C;
    end

    % now when ussing RxFFE floating taps need to tag stems correctly and
    % make sure DFEfloating tap code does not get exectuted

    %
else
    Cmod=C;%just us the FFE taps, C, passed for filtering
end
%%
%% filter the pulse response with the solved FFE
%  (now option to avoid this and just return Cmod for speed up)
if return_V
    Vfiltered=FFE( Cmod , param.RxFFE_cmx,spui, V );
else
    Vfiltered=[];
end