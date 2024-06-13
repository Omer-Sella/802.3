function [ bmax floating_tap_locations] = floating_taps_1sttest( hisi,N_b,N_bf,N_bg,N_bmax, bmaxg, COOP )
% Richard Mellitz:  04/23/2019
% hisi is the isi 1 ui/sample
% N_b number of fixed dfe taps
% N_bf number of floating taps per group
% N_bg number of floating tap groups. 1 2 or 3 right now
% N_bmax number of ui for the max reach of the floating taps
% bmaxg limit for the floating taps
% COOP = 1 co-optimize banks , -0 sequenatial optmization
%
%
% function to remove isi or add noise above bmaxg
if ~exist('COOP','var'), COOP=0;end
if iscolumn(hisi); hisi=hisi.';end
hsis_in=hisi;
% find all the reduction group taken N_bf at a time
% we are looking for the group when when remove yield the miminim isi, h, power
best_sigma=inf;best_ig1=-1;best_ig2=-1;best_ig3=-1;
% add on switch and loop for each potential group
switch N_bg
    case 0
        bmax=0;
        return
    case 1
        end1=N_bmax-N_bf;
        end2=N_b+1;
        end3=N_b+1;
    case 2
        end1=N_bmax-N_bf;
        end2=N_bmax-N_bf;
        end3=N_b+1;
    case 3
        end1=N_bmax-N_bf;
        end2=N_bmax-N_bf;
        end3=N_bmax-N_bf;
end
if COOP
    for ig1= N_b+1:end1 % now remove the 2nd group
        hcap=  hrem(hisi,ig1,N_bf,bmaxg)  ;
        % loop for 2rd group
        for ig2= N_b+1: end2
            hcap2= hrem(hcap,ig2,N_bf,bmaxg)  ;
            if N_bg < 2; hcap2 =hcap; end
            for ig3= N_b+1: end3
                hcap3= hrem(hcap2,ig3,N_bf,bmaxg)  ;
                if N_bg < 3 ; hcap3=hcap2 ; end
                sigma=norm( hcap3 );
                if sigma < best_sigma
                    best_sigma=sigma;
                    best_ig1=ig1;
                    best_ig2=ig2;
                    best_ig3=ig3;
                    best_hcap=hcap3;
                end
            end
        end
    end
else % sequentail
    for ig1= N_b+1:end1 % now remove the 1st group
        hcap=  hrem(hisi,ig1,N_bf,bmaxg)  ;
        sigma=norm( hcap );
        if sigma < best_sigma
            best_sigma=sigma;
            best_ig1=ig1;
            best_hcap=hcap;
        end
    end
    % loop for 2rd group
    hisi=best_hcap;
    for ig2= N_b+1: end2
        hcap= hrem(hisi,ig2,N_bf,bmaxg)  ;
        sigma=norm( hcap );
        if sigma < best_sigma
            best_sigma=sigma;
            best_ig2=ig2;
            best_hcap=hcap;
        end
    end
    hisi=best_hcap;
    % loop for 3rd group
    for ig3= N_b+1: end3
        hcap= hrem(hisi, ig3,N_bf,bmaxg)  ;
        sigma=norm( hcap );
        if sigma < best_sigma
            best_sigma=sigma;
            best_ig3=ig3;
            best_hcap=hcap;
        end
    end
    
end
bmax(N_b+1:N_bmax)=zeros(1,N_bmax-N_b);
switch N_bg
    case 1
        bmax(best_ig1:best_ig1+N_bf-1)=ones(1,N_bf)*bmaxg;
        floating_tap_locations= [best_ig1:best_ig1+N_bf-1];
    case 2
        bmax(best_ig1:best_ig1+N_bf-1)=ones(1,N_bf)*bmaxg;
        bmax(best_ig2:best_ig2+N_bf-1)=ones(1,N_bf)*bmaxg;
        floating_tap_locations= [best_ig1:best_ig1+N_bf-1 best_ig2:best_ig2+N_bf-1 ];
    case 3
        bmax(best_ig1:best_ig1+N_bf-1)=ones(1,N_bf)*bmaxg;
        bmax(best_ig2:best_ig2+N_bf-1)=ones(1,N_bf)*bmaxg;
        bmax(best_ig3:best_ig3+N_bf-1)=ones(1,N_bf)*bmaxg;
        floating_tap_locations= [best_ig1:best_ig1+N_bf-1 best_ig2:best_ig2+N_bf-1 best_ig3:best_ig3+N_bf-1 ];
end
floating_tap_locations=sort(floating_tap_locations);
if 0 % for code debug
    close force all
    stem(best_hcap,'disp','hcap')
    hold on
    stem(bmax,'-k','disp','bmax')
    stem(hisi,'disp','hisi')
    hold off
end
