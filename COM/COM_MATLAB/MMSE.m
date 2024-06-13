function MMSE_results = MMSE(PSD_results,sbr, cursor_i ,param, OP ) ;
if 1
    num_ui=param.num_ui_RXFF_noise;
    M=param.samples_per_ui;
    L=param.levels;
    sigma_X2=(L^2-1)/(3*(L-1)^2);
    fb=param.fb;
    R_LM=param.R_LM;
end
h=sbr(mod(cursor_i,M)+1:end-mod(cursor_i,M)); % align to sample point
h=reshape(h,1,[]); % make row vectors
h=[ h(1:floor(length(h)/M)*M) ];
h= [h zeros(1,num_ui*M-length(h)) ];
h=h(1:M:end);% resample
N=length(h);
dw=param.RxFFE_cmx ; % equalizer precuror tapsindx(1:N)=(1:N)-5-1;
dh=(cursor_i-mod(cursor_i,M))/M ; % precuror taps in h
if param.N_bg == 0
    Nw= param.RxFFE_cmx+1+param.RxFFE_cpx; %  total number of equalizer taps
    bmax=param.bmax;
    bmin=param.bmin ;
    wmax= [   ones(1,param.RxFFE_cmx-1)*param.ffe_tapn_max   param.ffe_pre_tap1_max 1.0  param.ffe_post_tap1_max     ones(1,param.RxFFE_cpx-1)*param.ffe_tapn_max  ];
    wmin= [  -ones(1,param.RxFFE_cmx-1)*param.ffe_tapn_max  -param.ffe_pre_tap1_max 1.0  -param.ffe_post_tap1_max   -ones(1,param.RxFFE_cpx-1)*param.ffe_tapn_max  ];
    idx=[];
else
    Nfloating_taps=param.N_bf*param.N_bg;
    Nmax=param.N_bmax;
    Nfix=param.RxFFE_cmx+1+param.RxFFE_cpx;
    Ng=param.N_bg;
    Nf=param.N_bf;
    Nw= dw+Nmax+1;%  total span of equalizer taps including floating taps
    Nwft=param.RxFFE_cmx+1+param.RxFFE_cpx+Nfloating_taps;%  total number of equalizer taps including floating taps
    hisi=h(dh+2:((dh-dw)+Nw));
    [idx]=findbankloc( hisi ,param.N_tail_start,param.N_bmax,Nf,inf,param.bmaxg,Ng  ); % using maximum power in hisi
    idx=sort(idx);
    bmax=param.bmax;
    bmin=param.bmin ;
    wmax= [   ones(1,param.RxFFE_cmx-1)*param.ffe_tapn_max   param.ffe_pre_tap1_max 1.0  param.ffe_post_tap1_max   ones(1,param.RxFFE_cpx-1)*param.ffe_tapn_max  ones(1,Nfloating_taps)*param.bmaxg ];
    wmin= [  -ones(1,param.RxFFE_cmx-1)*param.ffe_tapn_max  -param.ffe_pre_tap1_max 1.0  -param.ffe_post_tap1_max -ones(1,param.RxFFE_cpx-1)*param.ffe_tapn_max -ones(1,Nfloating_taps)*param.bmaxg  ];
    if 0
        figure
        set(gcf, 'tag', 'COM');%movegui(gcf,'south');
        hindx=1:length(hisi);
        stem(hindx,hisi)
        hold on
        stem(idx,hisi(idx))
    end

end
Nb=param.ndfe; % DFE taps
d=dw+dh; % used for index in algorithms
indx(1:N)=(1:N)-dh-1;

S_n=PSD_results.S_n; % total agregate noise PSD

Rn=ifft(S_n)*fb;
%% HH and R
Rnn=toeplitz(Rn(1:Nw),Rn(1:Nw));
hc1=[ h  zeros(1,Nw-1) ];
hr1=[ h(1) zeros(1,Nw-1)];
H=toeplitz(hc1,hr1);
if param.N_bg ~= 0
    H=H( :,[1:Nfix idx]);
end
HH= H'*H;
if param.N_bg ~= 0
 Rnn=Rnn( [1:Nfix idx-param.N_tail_start+1+Nfix],[1:Nfix idx-param.N_tail_start+1+Nfix]);
end
R=HH+Rnn/sigma_X2;
%% hb and h0
Hb= H(d+2:d+Nb+1,:);
h0=H(d+1,:);
% display(floor(h0));

%% Ib and zb (slide 10)
ib=eye(Nb);
zb=zeros(1,Nb);
wbl= [ R  -Hb' -h0';...
    -Hb  ib  zb'; ...
    h0  zb   0]\[h0'; zb' ;1];

%% re-adjust Nw to number of used taps
if param.N_bg ~= 0
    Nw=Nwft;
end
%% check equalized pulse
w=wbl(1:Nw);
b=wbl(Nw+1:length(wbl)-1); % dfe taps before limits are applied

%% apply blim (slide 11)  <---- need help here How do I get to RxFFE tap coefficents, C?


blim = min(bmax(:), max(bmin(:), b));
if (Nb > 0) && ~isequal(b, blim)
    wl = [R, -h0'; h0, 0]\[h0'+Hb'*blim; 1];
    w = wl(1:Nw);
end

wlim = min(wmax(:)*w(1+dw), max(wmin(:)*w(1+dw), w));
if ~isequal(w, wlim)
    wlim = wlim/(h0*wlim); % Ensure the equalized pulse amplitude is 1.
    if Nb > 0
        b = Hb*wlim; % Update the feedback coefficients.
        blim = min(bmax(:), max(bmin(:), b));
    end
    wl = [R, -h0'; h0, 0]\[h0'+Hb'*blim; 1];
    wl = wl(1:Nw);
    w = min(wmax(:)*wl(1+dw), max(wmin(:)*wl(1+dw), wl));
end

w=w(1:Nw) ;
% w=wl(1:Nw); % RxFFE taps before limits are applied
Craw=w/w(dw+1); % returned Rx FFE taps
% re-align Cmod to floating tap locations
if param.N_bg ~= 0
    C=Craw;
    C(Nfix+1:Nmax)=0;
    C(idx-param.N_tail_start+1+Nfix)=Craw(Nfix+(1:Nfloating_taps));
else
    C=Craw;
end
if 0
    figure
    set(gcf, 'tag', 'COM');movegui(gcf,'south');
     subplot(3,1,1)
     stem(-dw:Nmax-dw-1,C,'DisplayName','Returned RxFFE taps')
     title('solved Rxffe C')
    legend show
    subplot(3,1,2)
    resp=conv(C,h);
    resp=resp(dw+1:N+dw);
    stem(indx(1:length(resp)),resp,'disp', 'Rxffe filtered h with C')
    title('convolved response with C')
    legend show
    subplot(3,1,3)
    resp=conv(w,h);
    resp=resp(dw+1:N+dw);
    stem(indx(1:length(resp)),resp,'disp', 'Rxffe filtered h with w(ts) not normalized')
    title('convolved response with  w')
    legend show
end

MMSE_results.floating_tap_locations=idx;
MMSE_results.C=C;
MMSE_results.sigma_e=(w'*R*w+1+b'*b-2*w'*h0'-2*w'*Hb'*b);
MMSE_results.FOM=20*log10((R_LM/(L-1)/MMSE_results.sigma_e));

