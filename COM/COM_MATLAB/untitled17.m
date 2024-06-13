%pdf_full is self ISI pdf for each sample point
%h_j_full is the v/t calculation for each sample point
%the center vector for each should be identical to the standard COM variables
[pdf_full_1,h_j_full,A_s_vec] = get_pdf_full(chdata(1), delta_y, fom_result.t_s, param, OP,pdf_range) ;



if isempty(pdf_range)
    pdf_range=1:samp_UI;
else
    pdf_range=min(pdf_range):max(pdf_range);
end

%Test doing Level PDFs
Levels = 2*(0:param.levels-1)/(param.levels-1)-1;
%A_s_vec=A_s_vec*(param.levels-1)/param.R_LM;
A_s_vec=A_s_vec*(param.levels-1);

%add signal vector into pdf
for n=1:param.levels
    pdf_full{n}=pdf_full_1;
    for j=pdf_range
        pdf_full{n}(j).x=pdf_full{n}(j).x+A_s_vec(j)*Levels(n);
        pdf_full{n}(j).Min=pdf_full{n}(j).x(1)/pdf_full{n}(j).BinSize;
    end
end


% figure;
% hold on;
%This loop builds the same PDF/CDF structures from regular COM, but it is
%computed for every sample point
for n=1:param.levels
    for j=pdf_range
        sigma_G_full(j)  = norm([param.sigma_RJ*param.sigma_X*norm(h_j_full(:,j)), Struct_Noise.sigma_N, Struct_Noise.sigma_TX]);
        gaussian_noise_pdf_full(j) = normal_dist(sigma_G_full(j), Struct_Noise.ber_q, delta_y);
        gaussian_noise_pdf_full(j) = conv_fct(gaussian_noise_pdf_full(j), Struct_Noise.ne_noise_pdf);
        p_DD_full(j) = get_pdf_from_sampled_signal(param.A_DD*h_j_full(:,j), param.levels, delta_y);
        noise_pdf_full(j)=conv_fct(gaussian_noise_pdf_full(j), p_DD_full(j));
        isi_and_xtalk_pdf_full(j) = conv_fct_MeanNotZero(pdf_full{n}(j), Struct_Noise.cci_pdf);
        % change from adam gregory to include crosstalk
        %     combined_interference_and_noise_pdf_full = conv_fct(pdf_full(j), noise_pdf_full(j));
        combined_interference_and_noise_pdf_full{n}(j) = conv_fct_MeanNotZero(isi_and_xtalk_pdf_full(j), noise_pdf_full(j));
        
        %PDF to CDF
        combined_interference_and_noise_cdf_full{n}(j)=pdf_to_cdf(combined_interference_and_noise_pdf_full{n}(j));
        
    end
end
%hold off;


%For the given BER, find the top & bottom voltage level in the CDF
for n=1:param.levels
    A_ni_bottom{n}=zeros(1,samp_UI);
    A_ni_top{n}=zeros(1,samp_UI);
    for j=pdf_range
        [A_ni_top{n}(j),A_ni_bottom{n}(j)]=cdf_to_ber_contour(combined_interference_and_noise_cdf_full{n}(j),param.specBER);
    end
end
%plot(1:samp_UI,cursor_vector-A_ni_top,1:samp_UI,-cursor_vector+A_ni_bottom)

for n=1:param.levels-1
    eye_contour{n}(:,1)=A_ni_top{n+1};
    eye_contour{n}(:,2)=A_ni_bottom{n};
end


for n=1:param.levels-1
    %eye_contour holds the top eye in the 1st column & bottom eye in the 2nd column
    %define vref as middle of top eye height and bottom eye height.  Now
    %that all eyes are created, vref is non-zero except for middle eye
    EH_top=eye_contour{n}(half_UI,1);
    EH_bot=eye_contour{n}(half_UI,2);
    EH=EH_top-EH_bot;
    vref=EH_top/2+EH_bot/2;
    %This function finds left/right eye width by finding the vref crossings of
    %the top and bottom eye contours
    [Left_EW(n),Right_EW(n)]=find_eye_width(eye_contour{n},half_UI,samp_UI,vref);
end

%For reporting to .csv, need eye contour to be a matrix instead of cell
eye_contour_tmp=eye_contour;
eye_contour=[];
for n=1:param.levels-1
    eye_contour(:,(n-1)*2+1:n*2)=eye_contour_tmp{n};
end


%Find VEC eye height
out_VT=[];
out_VB=[];
if param.T_O ~=0
    
    switch lower(OP.Histogram_Window_Weight)
        case {'gaussian' 'norm' 'normal' 'guassian'}
            %build a gaussian window of weights that are multiplied by each pdf in the T_O range
            %Sigma = T_O/QL.  Default QL=2.5.  This gives a nice descent to near 0 at the edge of the window
            QL_sigma=T_O/param.QL;
            weights=exp(-1/2 * ([-T_O:T_O]/QL_sigma).^2);
        case 'triangle'
            %triangle window. linear slope from 0 to 1 and back down to 0
            %for the weights
            t_slope=1/(T_O);
            weights=[0:t_slope:1 1-t_slope:-t_slope:0];
        case 'rectangle'
            %default = rectangle.  all weights = 1
            weights(1:2*T_O+1)=1;
        case 'dual_rayleigh'
            QL_sigma=T_O/param.QL;
            X=-T_O:T_O;
            weights=(X+T_O)/QL_sigma^2.*exp(-1/2 * ((X+T_O)/QL_sigma).^2)...
                -(X-T_O)/QL_sigma^2.*exp(-1/2 * ((X-T_O)/QL_sigma).^2);
            weights=weights/max(weights);
        otherwise
            error('%s not recognized for Histogram_Window_Weight',OP.Histogram_Window_Weight)
    end

    for n=1:param.levels
        out_pdf{n}=[];
        for j=start_sample:end_sample
            target_pdf=combined_interference_and_noise_pdf_full{n}(j);
            target_pdf.y=target_pdf.y*weights(j-start_sample+1);
            if isempty(out_pdf{n})
                out_pdf{n}=target_pdf;
            else
                out_pdf{n} = combine_pdf_same_voltage_axis(out_pdf{n}, target_pdf);
            end
        end
        out_pdf{n}.y=out_pdf{n}.y/sum(out_pdf{n}.y);
    end
    
    for n=1:param.levels
        out_cdf{n}=pdf_to_cdf(out_pdf{n});
    end
    
    for n=1:param.levels
        [A_ni_top_O(n),A_ni_bottom_O(n)]=cdf_to_ber_contour(out_cdf{n},param.specBER);
    end
    
    for n=1:param.levels-1
        OUT_VT_L(n,1)=A_ni_top_O(n+1);
        OUT_VT_L(n,2)=A_ni_bottom_O(n);
    end
    
    %Report the top/bottom eye height of the worst eye
    EH_VT=OUT_VT_L(:,1)-OUT_VT_L(:,2);
    [mineh,min_idx]=min(EH_VT);
    out_VT=OUT_VT_L(min_idx,1);
    out_VB=OUT_VT_L(min_idx,2);
    
%     CDF_Mean=mean(A_s_vec(start_sample:end_sample));
%     out_VT=2*CDF_Mean-A_ni_top_O;
%     out_VB=-1*A_ni_bottom_O;
    
    if debug_plot
        figure;
        hold on;
        for j=start_sample:end_sample
            plot(combined_interference_and_noise_pdf_full(j).x,combined_interference_and_noise_pdf_full(j).y);
        end
        plot(out_pdf.x,out_pdf.y,'color','k','LineWidth',2);
        hold off;
    end
end
