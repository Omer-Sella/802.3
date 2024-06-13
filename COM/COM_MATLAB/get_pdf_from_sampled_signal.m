function [ pdf ] = get_pdf_from_sampled_signal( input_vector, L, BinSize ,FAST_NOISE_CONV)
% Create PDF from interference vector using successive delta-set convolutions.
%   input_vector = list of values of samples
%   return
%   pdf.x
%   pdf.y
%   pdf.vec
%   pdf.bin
if ~exist('FAST_NOISE_CONV','var')
    FAST_NOISE_CONV=0;
end
if max(input_vector) > BinSize
    input_vector=input_vector(abs(input_vector)>BinSize);
end
% for i = 1:length(input_vector)
%    if abs(input_vector(i)) < BinSize , input_vector(i)=0; end
%end

input_vector(abs(input_vector)<BinSize) = 0;
b=sign(input_vector);
[input_vector,index]=sort(abs(input_vector),'descend');
input_vector=input_vector.*b(index);
if FAST_NOISE_CONV
    sig_res=norm(input_vector(find(abs(input_vector)<.001,1)+1:end));
    res_pdf= normal_dist(sig_res,5,BinSize);
    input_vector=input_vector(1:find(abs(input_vector)<.001,1));
end
%% Equation 93A-39 %%
values = 2*(0:L-1)/(L-1)-1;
prob = ones(1,L)/L;

%% Initialize pdf to delta at 0
pdf=d_cpdf(BinSize, 0, 1);
empty_pdf=pdf;
for k = 1:length(input_vector)
    %     pdfn=d_cpdf(BinSize, abs(input_vector(k))*values, prob);
    pdfn=Init_PDF_Fast(empty_pdf, abs(input_vector(k))*values, prob);
    pdf=conv_fct(pdf, pdfn);
end
if FAST_NOISE_CONV
%     pdf=conv_fct(pdf,res_pdf);
    pdf=conv_fct_TEST(pdf,res_pdf);
end

function [pdf,h_j_full,A_s_vec]=get_pdf_full(chdata, delta_y, t_s, param, OP,pdf_range)
t_s_orig=t_s;
%SBR=chdata.eq_pulse_response(:)'; % row vector SBR(1:t_s-14)=0;SBR(t_s+15:end)=0; for debug (AJG edit)
type=chdata.type;

pulse_orig=chdata.eq_pulse_response(:)';
%build arbitrary time axis with step size = 1/samples per ui
old_time=[0:length(pulse_orig)-1]/param.samples_per_ui;
%force t_s at time =0 (makes the other things below easy)
original_sample_time=old_time(t_s_orig);
old_time=old_time-original_sample_time;
%build new time axis that forces time=0 to be in the axis
%unless the new/old samples per UI are integer ratios, time 0 will not be
%there by default
samp_UI=param.samples_for_C2M;
new_timea=[0:-1/samp_UI:min(old_time)];
new_timeb=[0:1/samp_UI:max(old_time)];
new_time=[fliplr(new_timea) new_timeb(2:end)];
SBR=interp1(old_time,pulse_orig,new_time);
%new sample time is simply the point where new_time = 0
[tmp,t_s]=min(abs(new_time));

residual_response = SBR;

half_UI=get_center_of_UI(samp_UI);

if isequal(type, 'THRU')
    % for thru pulse response:
    % remove the cursor and the DFE postcursors (up to their limit), since
    % we only care about the residuals.
    
    %AJG021820
    if ~param.Floating_DFE
        ideal_cancelled_cursors = SBR(t_s+samp_UI*(1:param.ndfe));
    else
        ideal_cancelled_cursors = SBR(t_s+samp_UI*(1:param.N_bmax));
    end
    if param.dfe_delta ~= 0
        ideal_cancelled_cursors_q=floor(abs( ideal_cancelled_cursors/(residual_response(t_s).*param.dfe_delta) )).*residual_response(t_s)*param.dfe_delta.*sign(ideal_cancelled_cursors);
    else
        ideal_cancelled_cursors_q=ideal_cancelled_cursors;
    end
    
    if ~param.Floating_DFE
        bmax_vec=residual_response(t_s)*[param.bmax];
        bmin_vec=residual_response(t_s)*[param.bmin];
    else
        bmax_vec=residual_response(t_s)*[param.use_bmax];
        bmin_vec=residual_response(t_s)*[param.use_bmin];
    end
    effective_cancelled_cursors=dfe_clipper(ideal_cancelled_cursors_q,bmax_vec,bmin_vec);
    
    
    effective_cancellation_samples = kron(effective_cancelled_cursors, ones(1, samp_UI));
    dfetaps=effective_cancelled_cursors/SBR(t_s);
    
    % Apply a constant DFE coefficient 1/2 UI before and after each postcursor. Not
    % really needed for COM, but helps debugging. May be factored out in future revisions.
    
    %avoid dividing samp_UI by 2 in case it is not even
    start_cancel=t_s-half_UI+1+samp_UI;
    %AJG021820
    if ~param.Floating_DFE
        end_cancel=start_cancel+param.ndfe*samp_UI-1;
    else
        end_cancel=start_cancel+param.N_bmax*samp_UI-1;
    end
    residual_response(start_cancel:end_cancel) = ...
        residual_response(start_cancel:end_cancel) - effective_cancellation_samples;
    %else
    % for crosstalk pulse responses, nothing is cancelled, and all phases
    % are equally important.
    
    %remove entire cursor UI
    uiv_start=start_cancel-samp_UI;
    uiv_end=uiv_start+samp_UI-1;
    A_s_vec = param.R_LM*SBR(uiv_start:uiv_end)/(param.levels-1);
    residual_response(uiv_start:uiv_end)=0;
end

nui=round(length(residual_response)/samp_UI);


vs=transpose(reshape(residual_response(samp_UI+1:samp_UI*(nui-1)),samp_UI,nui-2));
%added vs_raw in order to calculate h_j.  vs_raw uses the pulse
%response without DFE included.  (Can't include DFE for jitter calc)
vs_raw=transpose(reshape(SBR(samp_UI+1:samp_UI*(nui-1)),samp_UI,nui-2));

% if OP.DISPLAY_WINDOW,
%     hwaitbar=waitbar(0);
% end

% determine which pdf to use
if isequal(type, 'THRU')
    % one phase is interesting for thru
    phases = mod(t_s,samp_UI);
    if phases==0, phases = samp_UI; end
else
    phases=1:samp_UI;
end

mxV = zeros(size(phases));

%phases reveals the raw position in the UI window of the cursor.
%shift_amount is the amount to shift so that it aligns with half_UI
shift_amount=half_UI-phases;
%vs_shift puts the cursor at the center
vs_shift=circshift(vs,[0 shift_amount]);
L=size(vs_raw,1);
%allow partial UI computation through pdf_range
%if pdf_range is empty, do full UI
if isempty(pdf_range)
    pdf_range=1:samp_UI;
else
    pdf_range=min(pdf_range):max(pdf_range);
end
h_j_full=zeros(L,samp_UI);
for k=pdf_range
    pdf(k)=get_pdf_from_sampled_signal(vs_shift(:,k), param.levels, delta_y); %#ok<AGROW>
    %mxV(k)=sqrt(sum( pdf_samples(k).x.^2.*pdf_samples(k).y)); % standard deviation of PDF
    %progress = k/length(phases);
    %if OP.DISPLAY_WINDOW, waitbar(progress, hwaitbar, ['processing COM pdf ' chdata.base ] ); figure(hwaitbar); drawnow; end
    
    %build the circshift of h_j_full into the loop to support a reduced
    %range of sampling points. circshift at the end only works if doing the
    %full range of sampling points.  And shifting before the loop will
    %yield the wrong answer at the edges of the UI
    hk=k-shift_amount;
    if hk<1
        hk=hk+samp_UI;
    elseif hk>samp_UI
        hk=hk-samp_UI;
    end
    if hk==1
        %when hk=1, the early UI is the last column
        h_j_full(1:L-1,k)=(vs_raw(2:end,hk+1)-vs_raw(1:end-1,samp_UI))/2*samp_UI;
    elseif hk==samp_UI
        %when hk=samp_UI, the late UI is the first column
        h_j_full(1:L-1,k)=(vs_raw(2:end,1)-vs_raw(1:end-1,hk-1))/2*samp_UI;
    else
        %for all other cases, do the normal late=+1, early = -1
        h_j_full(1:L,k)=(vs_raw(:,hk+1)-vs_raw(:,hk-1))/2*samp_UI;
    end
end