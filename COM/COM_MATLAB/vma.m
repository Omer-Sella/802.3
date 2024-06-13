function results= vma(PR, M)
% PR=sbr.Data;
% M=32;
% PR is the pulse response
% M is samples per UI
[ seq, syms, syms_nrz ] = PRBS13Q( );
% seq uses [ -1 -1/3 1/3 1] & syms uses [ 0 1 2 3 ]
symbols=seq; 
imaxPR=find(PR==max(PR),1,'first'); % find index for peak
% start end symbols index for 7 3's and 6 0's
indx_S3x7_start=M*(strfind(syms,ones(1,7)*3)-1)+imaxPR;
indx_S3x7_end=M*(strfind(syms,ones(1,7)*3)+5)+imaxPR;
indx_S0x6_start=M*(strfind(syms,ones(1,6)*0))-1+imaxPR;
indx_S0x6_end=M*(strfind(syms,ones(1,6)*0)+4)+imaxPR;
% superposition code
shifting_vector=kron(symbols,[ 1  zeros(1,M-1) ]) ;
Bit_stream_response=filter(PR,1, shifting_vector);
% find center of 3's and 0's
icent3=floor((indx_S3x7_end-indx_S3x7_start)/2 + indx_S3x7_start);
icent0=floor((indx_S0x6_end-indx_S0x6_start)/2 + indx_S0x6_start);
% plot(Bit_stream_response(indx_S3x7_start:indx_S3x7_end))
% hold on
% plot(Bit_stream_response(indx_S0x6_start:indx_S0x6_end))
P_3 = mean( Bit_stream_response((icent3-M):(icent3+M) ) );
P_0 = mean( Bit_stream_response((icent0-M):(icent0+M) ) );
VMA= P_3 - P_0;
results.P_3=P_3;
results.P_0=P_0;
results.VMA=VMA;
function line_intersection=vref_intersect(eye_contour,x_in,vref)

%slope of the 2 sample points around vref crossing
m1=(eye_contour(x_in,1)-eye_contour(x_in-1,1));
%x-intercept for the line
b1=eye_contour(x_in,1)-m1*x_in;
% drawing a horizontal line through vref so slope = 0
m2=0;
%special case for horizontal line, b=y
b2=vref;
%the x-value of line intersection = (b2-b1)/(m1-m2)
%sinc m2 is always 0 and b2 is always vref, this could be stated as (vref-b1)/m1
%And usually vref is 0, so it further reduces to -b1/m1
line_intersection=(b2-b1)/(m1-m2);

