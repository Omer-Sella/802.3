function [sigma_XT, sigma_FEXT, sigma_NEXT] = get_xtlk_noise( upsampled_txffe, type, param, chdata, phase_memory,C )
% Modified not to double count crosstalk: John Ewen 13/12/2018
% function sigma_XT = get_xtlk_noise( upsampled_txffe, type, param, chdata,ctle_indx,clow_indx, C,cursor_i)
index_f2=find(chdata(1).faxis(:)>param.fb,1,'first');
if isempty(index_f2), index_f2=length(chdata(1).faxis);end
f=chdata(1).faxis;
temp_angle=(param.samples_per_ui*param.sample_dt)*pi.*chdata(1).faxis;
if(f(1)==0)
    temp_angle(1)=1e-20;% we don't want to divide by zero
end
PWF_tx=ones(1,length(f));
if max(upsampled_txffe) > 0
    PWF_tx=zeros(1,length(f));
    [mcur,icur] = max(upsampled_txffe);
    if exist('phase_memory','var') && ~isempty(phase_memory)
        pre_calc=1;
    else
        pre_calc=0;
    end
    for ii=1:length(upsampled_txffe)
        if upsampled_txffe(ii)==0
            %speed up:  skip cases when txffe=0
            continue;
        end
        %         PWF_tx=upsampled_txffe(ii).*exp(-1j*2*pi*(ii-icur).*f/param.fb)+PWF_tx;
        % Adee Ran 2020-06-03 remove create large 2D matrix when 1D is all
        % that is needed
        if ii==icur
            %speed up:  ii-icur=0, so just scalar addition and avoid exp calc
            PWF_tx = PWF_tx + upsampled_txffe(ii);
        else
            if pre_calc
                %speed up:  avoid vector exp calculation by externally pre-calculating it
                term_ii = upsampled_txffe(ii).*phase_memory(:,ii);
            else
                term_ii = upsampled_txffe(ii).*exp(-1j*2*pi*(ii-icur).*f/param.fb);
            end
            %bug fix:  use transpose instead of ' to avoid taking complex conjugate
            PWF_tx = PWF_tx + transpose(term_ii(:));
        end
        % /Adee
    end
end
PWF_rx=ones(1,length(f));
if exist('C','var')
    PWF_rx=zeros(1,length(f));
    for ii=-param.RxFFE_cmx:param.RxFFE_cpx
        if C(ii+param.RxFFE_cmx+1)==0
            %speed up:  skip cases when rxffe=0
            continue;
        end
        if ii+1==0
            %speed up:  ii+1=0, so just scalar addition and avoid exp calc
            PWF_rx = PWF_rx + C(ii+param.RxFFE_cmx+1);
        else
            if pre_calc
                %speed up:  avoid vector exp calculation by externally pre-calculating it
                %The latter columns of phase_memory hold RXFFE shift vectors
                term_ii=C(ii+param.RxFFE_cmx+1).*phase_memory(:,ii+param.RxFFE_cmx+1+length(upsampled_txffe));
                term_ii=transpose(term_ii);
            else
                term_ii=C(ii+param.RxFFE_cmx+1).*exp(-1j*2*pi*(ii+1).*f/param.fb);
            end
            PWF_rx=PWF_rx+term_ii;
        end
        %PWF_rx=C(ii+param.RxFFE_cmx+1).*exp(-1j*2*pi*(ii+1).*f/param.fb)+PWF_rx;
    end
end
MDFEXT=0;MDNEXT=0;MDNEXT_ICN=0;MDFEXT_ICN=0;
for ii=2:size(chdata,2)
    SINC = sin(temp_angle)./temp_angle;
    PWF_data=SINC.^2;
    PWF=PWF_data.*PWF_rx; % power weight function
    PWFnext=abs(PWF);
    PWF=PWF_data.*PWF_rx.*PWF_tx; % power weight function
    PWFfext=abs(PWF);
    if isequal(chdata(ii).type, 'FEXT')
        MDFEXT=sqrt(abs(chdata(ii).sdd21ctf).^2+MDFEXT.^2); % power sum xtk
        MDFEXT_ICN=sqrt(2*chdata(ii).delta_f/param.f2*sum( chdata(ii).A^2*PWFfext(1:index_f2).*abs(MDFEXT(1:index_f2)).^2)); %eq 46
    elseif isequal(chdata(ii).type, 'NEXT')
        MDNEXT=sqrt(abs(chdata(ii).sdd21ctf).^2+MDNEXT.^2); % power sum xtk
        MDNEXT_ICN=sqrt(2*chdata(ii).delta_f/param.f2*sum( chdata(ii).A^2*PWFnext(1:index_f2).*abs(MDNEXT(1:index_f2)).^2)); %eq 47
    end
end
if nargout == 1 && isequal(type,'NEXT')
    sigma_XT = MDNEXT_ICN*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
elseif nargout == 1 && isequal(type,'FEXT')
    sigma_XT = MDFEXT_ICN*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
elseif nargout == 3
    sigma_XT=norm([ MDNEXT_ICN MDFEXT_ICN ])*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
    sigma_NEXT = MDNEXT_ICN*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
    sigma_FEXT = MDFEXT_ICN*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
end
