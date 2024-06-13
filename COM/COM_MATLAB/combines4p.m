function [ s11out, s12out, s21out, s22out ] = combines4p( s11in1, s12in1, s21in1, s22in1, s11in2, s12in2, s21in2, s22in2)

%original method:
% s1=zeros(2,2,length(s11in1)); s2=s1; t3=s1;
% for i=1:length(s11in1)
%     s1(:,:,i)=[s11in1(i) s12in1(i); s21in1(i) s22in1(i) ];
%     s2(:,:,i)=[s11in2(i) s12in2(i); s21in2(i) s22in2(i) ];
% end
% t1=stot(s1);
% t2=stot(s2);
% for i=1:length(s11in1)
%     t3(:,:,i)=t1(:,:,i)*t2(:,:,i);
% end
% s3=ttos(t3);
% s11out=s3(1,1,:);
% s11out=transpose(s11out(:));
% s12out=s3(1,2,:);
% s12out=transpose(s12out(:));
% s21out=s3(2,1,:);
% s21out=transpose(s21out(:));
% s22out=s3(2,2,:);
% s22out=transpose(s22out(:));


%vectorized method:
s1(1,1,:)=s11in1;
s1(1,2,:)=s12in1;
s1(2,1,:)=s21in1;
s1(2,2,:)=s22in1;
s2(1,1,:)=s11in2;
s2(1,2,:)=s12in2;
s2(2,1,:)=s21in2;
s2(2,2,:)=s22in2;


N = (1-s1(2,2,:).*s2(1,1,:)) ;
s11out = s1(1,1,:)+(s1(1,2,:).*s1(2,1,:).*s2(1,1,:))./N ;
s12out = s1(1,2,:).*s2(1,2,:)./N ;
s21out = s2(2,1,:).*s1(2,1,:)./N;
s22out = s2(2,2,:)+(s2(1,2,:).*s2(2,1,:).*s1(2,2,:))./N ;

s11out=transpose(squeeze(s11out));
s12out=transpose(squeeze(s12out));
s21out=transpose(squeeze(s21out));
s22out=transpose(squeeze(s22out));
function p=conv_fct(p1, p2)
if p1.BinSize ~= p2.BinSize
    error('bin size must be equal')
end

p=p1;
%p.BinSize=p1.BinSize;
%p.Min=p1.Min+p2.Min;
p.Min=round(p1.Min+p2.Min);	% modified by Yasuo Hidaka, 9/4/2016
p.y=conv2(p1.y, p2.y);
%p.x =p.Min*p.BinSize:p.BinSize:-p.Min*p.BinSize;
%p.x =(p.Min:-p.Min)*p.BinSize;	% modified by Yasuo Hidaka, 9/4/2016
pMax=p.Min+length(p.y)-1;
p.x =(p.Min*p.BinSize:p.BinSize:pMax*p.BinSize);
