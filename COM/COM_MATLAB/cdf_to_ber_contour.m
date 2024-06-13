function [noise_bottom,noise_top]=cdf_to_ber_contour(cdf,specBER)


%For the given BER, find the top & bottom voltage level in the CDF

%for the top, just find the first index at the spec BER
nidx=find(cdf.y>specBER, 1, 'first');
noise_bottom = cdf.x(nidx);
%for top, flip the cdf.  need to operate on row vector for fliplr:  cdf.y(:)'
nidx=find(fliplr(cdf.y(:)')>specBER, 1, 'first');
%the true index without flipping
nidx=length(cdf.y)-nidx+1;
noise_top = cdf.x(nidx);
function p=comb_fct(p1, p2)
if p1.BinSize ~= p2.BinSize
    error('bin size must be equal')
end

p=p1;
p.BinSize=p1.BinSize;
%p.Min=p1.Min+p2.Min;
p.Min=min(p1.Min,p2.Min);
difsz=abs(p1.Min-p2.Min);
lp1=length(p1.y);
lp2=length(p2.y);
if p1.Min == p.Min
    p2.y(difsz+1:lp2+difsz)=p2.y;
    p2.y(1:difsz)=0;
    p2.y(lp2+difsz+1:lp1)=0;
elseif p2.Min == p.Min
    p1.y(difsz+1:lp1+difsz)=p1.y;
    p1.y(1:difsz)=0;
    p1.y(lp1+difsz+1:lp2)=0;
end
    % p.y=conv(p1.y, p2.y);
p.y=(p1.y+p2.y);
% p.y=p.y/sum(p.y);
% p.x =p.Min*p.BinSize:p.BinSize:-p.Min*p.BinSize;
p.x =(p.Min:-p.Min)*p.BinSize;	% modified by Yasuo Hidaka, 9/4/2016
