function [Left_EW,Right_EW,eye_contour,out_VT,out_VB]=COM_eye_width(chdata,delta_y,fom_result,param,OP,Struct_Noise,pdf_range_flag)


debug_plot=0;


samp_UI=param.samples_for_C2M;
half_UI=get_center_of_UI(samp_UI);
T_O=floor((param.T_O/1000)*param.samples_for_C2M);
start_sample=half_UI-T_O;
end_sample=half_UI+T_O;


%pdf_range is a placeholder to allow fractional UI calculation for optimize
%C2M speedup.  For regular COM_eye_width calls, pdf_range will be empty
%which will force pdf_range=1:samp_UI for stanard full UI calculation.
if pdf_range_flag
    pdf_range=[start_sample end_sample]; 
else
    pdf_range=[];
end