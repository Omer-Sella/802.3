function [Left_EW,Right_EW]=find_eye_width(eye_contour,half_UI,samples_per_UI,vref)

%Left eye Width (Top Eye)
left_top=eye_contour(half_UI:-1:1,1);
%vref_crossing is the first point less than vref (usually first point < 0)
vref_crossing=find(left_top<vref,1,'first');
if isempty(vref_crossing)
    %this case handles completely open eye
    L1=half_UI;
elseif vref_crossing==1
    %this case handles completely closed eye
    L1=0;
else
    %this case handles the normal eye
    %INT is a linear interpolation between the 2 points on either side of
    %vref to determine where the vref crossing occurred.  In systems with
    %a small number of samples_per_UI, interpolation improves accuracy over
    %just using the integer sample point
    INT=vref_intersect(eye_contour(1:half_UI,1),half_UI-vref_crossing+1+1,vref);
    L1=half_UI-INT;
end
%Left eye Width (Bottom Eye)
left_bot=eye_contour(half_UI:-1:1,2);
vref_crossing=find(left_bot>vref,1,'first');
if isempty(vref_crossing)
    L0=half_UI;
elseif vref_crossing==1
    L0=0;
else
    INT=vref_intersect(eye_contour(1:half_UI,2),half_UI-vref_crossing+1+1,vref);
    L0=half_UI-INT;
end
%Right eye Width (Top Eye)
right_top=eye_contour(half_UI:end,1);
vref_crossing=find(right_top<vref,1,'first');
if isempty(vref_crossing)
    R1=samples_per_UI-half_UI;
elseif vref_crossing==1
    R1=0;
else
    INT=vref_intersect(eye_contour(half_UI:end,1),vref_crossing,vref)+half_UI-1;
    R1=INT-half_UI;
end
%Right eye Width (Bottom Eye)
right_bot=eye_contour(half_UI:end,2);
vref_crossing=find(right_bot>vref,1,'first');
if isempty(vref_crossing)
    R0=samples_per_UI-half_UI;
elseif vref_crossing==1
    R0=0;
else
    INT=vref_intersect(eye_contour(half_UI:end,2),vref_crossing,vref)+half_UI-1;
    R0=INT-half_UI;
end

%L1 = top left eye width
%L0 = bottom left eye width
%Left eye width is the minimum
%R1 = top right eye width
%R0 = bottom right eye width
%Right eye width is the minimum
Left_EW=min([L1 L0]);
Right_EW=min([R1 R0]);