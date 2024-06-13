function  [chdata, param, SDDch,SDDp2p] = TD_FD_fillin(param, OP, chdata)
Over_sample=2;
num_files=length(chdata);
for i=1:num_files
    V=chdata(i).uneq_pulse_response;
    T=chdata(i).t;
    dt=T(2)-T(1);
    f=0:1/max(T):1/dt;
    chdata(i).faxis=f;
    f75=find(f >= param.fb*.75,1,'first');
    fnq=find(f >= param.fb*.5,1,'first');
    chdata(i).fmaxi = length(f);
    chdata(i).faxis = f;
    UI=param.ui; % unit interval
    M=param.samples_per_ui; % sample per UI
    N_v=param.N_v; % number of UI for Vf determination
    
    % filters
    H_bt=Bessel_Thomson_Filter(param,f,OP.Bessel_Thomson);
    H_bw=Butterworth_Filter(param,f,OP.Butterworth);
    H_ftr=H_bw.*H_bt;
    H_ftr=H_ftr(:);
    % fd of PR
    prr=sin(f*UI*pi)./(f*UI*pi); %nremoving sinc function to avoid using sig proc toolbox
    prr = prr(:);
    if f(1)==0, prr(1)=1; end %remove NaN
    fd=fft(V);
    fd=fd(1:floor(length(fd)/2)); % un process freq domain response
    
    %% get Vf
    shifting_vector=kron(ones(1,200) ,[ 1  zeros(1,M-1) ]) ;
    step_response=filter(V,1, shifting_vector);
    Vf=step_response(end);
    STEP=timeseries(step_response(1:length(shifting_vector)),T(1:length(shifting_vector)));
    %%
    % ILest=20.*log10(abs(fd(1:f75)/Vf/M/Over_sample))  -   20.*log10(abs(prr(1:f75))) - 20*log10(abs(squeeze(H_ftr(1:f75)))) ;  %removing db function to avoid using sig proc toolbox
    % figure
    % plot(f(1:f75),ILest)
    
    IL_conv=fd(1:f75)/Vf/M/Over_sample  ./  prr(1:f75) ./ H_ftr(1:f75) ;  %removing db function to avoid using sig proc toolbox
    % set same variables as get_s4p_files
    IL_fields={'sdd12_raw' 'sdd21_raw' 'sdd12_orig' 'sdd21_orig' 'sdd12' 'sdd21' 'sdd21p' 'sdd21f'};
    Zero_fields={'sdd22_raw' 'sdd11_raw' 'sdc12_raw' 'sdc21_raw' 'sdc22_raw' 'sdc11_raw' ...
        'sdd11_orig' 'sdd22_orig' 'sdd11' 'sdd22' 'sdc12' 'sdc21' 'sdc11' 'sdc22' 'sdc21p'};
    zero_vector=zeros(length(IL_conv),1);
    for j=1:length(IL_fields)
        chdata(i).(IL_fields{j})=IL_conv;
    end
    for j=1:length(Zero_fields)
        chdata(i).(Zero_fields{j})=zero_vector;
    end
    
    if i==1
        SDDch(:,1,2)=chdata.sdd12_raw;
        SDDch(:,2,1)=chdata.sdd21_raw;
        SDDch(:,1,1)=chdata.sdd11_raw;
        SDDch(:,2,2)=chdata.sdd22_raw;
        SDDp2p= zeros(length(IL_conv),1);
    end
    chdata(i).TX_RL=[];
    chdata(i).TDR11=[];
    chdata(i).PDTR11=[];
    chdata(i).TDR22=[];
    chdata(i).PDTR22=[];  
end