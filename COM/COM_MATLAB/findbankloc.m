function [idx]=findbankloc(hisi,idx_st,idx_en,tap_bk,curval,bmaxg,N_bg)
% [idx]=findbankloc(hisi,idx_st,idx_en,tap_bk)
% find the location of the DFE bank
% hisi: waveform with cursor values;
% idx_st: starting index;
% idx_en: ending index ;
% tap_bk: number of taps per bank;
% bmaxg: maximum coefficient;

hisi=hisi(:);
len=idx_en-idx_st+1;
h0=abs(hisi(idx_st:idx_en));
h1=max(0,h0-bmaxg*curval);

h0n=zeros(len-tap_bk+1,1);
h1n=h0n;

for ii=1:tap_bk
    h0tmp=h0(ii:ii+len-tap_bk);
    h0n=h0n+h0tmp.^2;
    h1tmp=h1(ii:ii+len-tap_bk);
    h1n=h1n+h1tmp.^2;
end

ndiff=h0n-h1n;



idx=zeros(1,tap_bk*N_bg);
ordered_set=(1:(N_bg-1)*tap_bk+1)';
set_next_bank=0;
%Loop through each group
for k=1:N_bg
    %Sort to choose the strongest
    [~,val_sort]=sort(ndiff,'descend');
    if k==1
        %shortcut:  Choose the first 1:N_bg*tap_bk taps if they are the strongest
        if isequal(sort(val_sort(ordered_set)),ordered_set)
            idx=1:N_bg*tap_bk;
            break;
        end
    end
    if set_next_bank>0
        %when a previous bank (goodV) was found, automatically set the bank without going through the search
        new_bank=set_next_bank:set_next_bank+tap_bk-1;
        idx(tap_bk*(k-1)+1:tap_bk*k)=new_bank;
        set_next_bank=0;
        ndiff(new_bank)=0;
        bad_start=new_bank(1)-tap_bk+1;
        bad_end=new_bank(1)-1;
        if bad_end<=0
            badV=[];
        elseif bad_start>0
            badV=bad_start:bad_end;
        else
            badV=1:bad_end;
        end
        if ~isempty(badV)
            ndiff(badV)=0;
        end
        continue;
    end
    %potential bank = the strongest tap group
    new_bank=val_sort(1):val_sort(1)+tap_bk-1;
    if k==N_bg
        %Last group:  just choose the strongest
        idx(tap_bk*(k-1)+1:tap_bk*k)=new_bank;
        break;
    end
    
    do_it_again=1;
    first_time=1;
    while do_it_again
        do_it_again=0;
        %note badV:  taps smaller and less than 1 group away
        bad_start=new_bank(1)-tap_bk+1;
        bad_end=new_bank(1)-1;
        if bad_end<=0
            badV=[];
        elseif bad_start>0
            badV=bad_start:bad_end;
        else
            badV=1:bad_end;
        end
        for j=length(badV):-1:1
            if any(badV(j)-idx==0)
                badV(j)=[];
            end
        end
        %note goodV:  the tap exactly 1 tap_bk smaller
        goodV=new_bank(1)-tap_bk;
        if ~isempty(badV)
            if ~first_time
                [~,val_sort]=sort(ndiff,'descend');
            end
            first_time=0;
            checkV=[badV new_bank];

            badV_pos=zeros(1,length(badV));
            for j=1:length(badV)
                badV_pos(j)=find(badV(j)==val_sort);
            end
            
            %loop through the sorted list to find the first tap outside the group and not a member of badV
            found_goodV=0;
            for ii=1:length(val_sort)
                if val_sort(ii)==goodV
                    found_goodV=1;
                    break;
                end
                if all(val_sort(ii)-checkV~=0)
                    break;
                end
            end
            
            if ~found_goodV && min(badV_pos)<ii
                %if goodV wasn't found and bad taps occur before non group members are found
                %throw out the strongest tap and take the next strongest
                do_it_again=1;
                ndiff(new_bank(1))=0;
                %speed up:  new max_val is always val_sort(2)
                new_bank=val_sort(2):val_sort(2)+tap_bk-1;
            end
            if found_goodV
                %if goodV was found, set the next bank to goodV
                set_next_bank=goodV;
            end
        end
        
    end
    %at the end, the floating taps are set to idx
    %and ndiff has illegal values set to zero
    ndiff(new_bank)=0;
    idx(tap_bk*(k-1)+1:tap_bk*k)=new_bank;
    if ~isempty(badV)
        ndiff(badV)=0;
    end
end


idx=idx+idx_st-1;