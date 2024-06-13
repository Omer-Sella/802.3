function pdf=get_pdf(chdata, delta_y, t_s, param, OP,ixphase)
SBR=chdata.eq_pulse_response(:)'; % row vector
type=chdata.type;
samp_UI=param.samples_per_ui;
residual_response = SBR;

if isequal(type, 'THRU')
    % for thru pulse response:
    % remove the cursor and the DFE postcursors (up to their limit), since
    % we only care about the residuals.
    
    if ~param.Floating_DFE
        ideal_cancelled_cursors = SBR(t_s+param.samples_per_ui*(0:param.ndfe));
    else
        ideal_cancelled_cursors = SBR(t_s+param.samples_per_ui*(0:param.N_bmax));
    end
    if param.dfe_delta ~= 0
        ideal_cancelled_cursors_q=floor(abs( ideal_cancelled_cursors/(residual_response(t_s).*param.dfe_delta) )).*residual_response(t_s)*param.dfe_delta.*sign(ideal_cancelled_cursors);
    else
        ideal_cancelled_cursors_q=ideal_cancelled_cursors;
    end
    
    %AJG021820
    if ~param.Floating_DFE
        bmax_vec=residual_response(t_s)*[1,param.bmax];
        bmin_vec=residual_response(t_s)*[1,param.bmin];
    else
        bmax_vec=residual_response(t_s)*[1,param.use_bmax];
        bmin_vec=residual_response(t_s)*[1,param.use_bmin];
    end
    effective_cancelled_cursors=dfe_clipper(ideal_cancelled_cursors_q,bmax_vec,bmin_vec);
    
    
    effective_cancellation_samples = kron(effective_cancelled_cursors, ones(1, param.samples_per_ui));
    dfetaps=effective_cancelled_cursors/SBR(t_s);
    
    % Apply a constant DFE coefficient 1/2 UI before and after each postcursor. Not
    % really needed for COM, but helps debugging. May be factored out in future revisions.
    start_cancel = t_s-param.samples_per_ui/2;
    if ~param.Floating_DFE
        end_cancel = t_s+(1/2+param.ndfe)*param.samples_per_ui - 1;
    else
        end_cancel = t_s+(1/2+param.N_bmax)*param.samples_per_ui - 1;
    end
    residual_response(start_cancel:end_cancel) = ...
        residual_response(start_cancel:end_cancel) - effective_cancellation_samples;
    %else
    % for crosstalk pulse responses, nothing is cancelled, and all phases
    % are equally important.
end

nui=round(length(residual_response)/param.samples_per_ui);

vs=zeros(nui-2, param.samples_per_ui);
for i=1:param.samples_per_ui
    vs(:,i)=residual_response(param.samples_per_ui*(1:nui-2)+i);
end

if OP.DISPLAY_WINDOW,
    hwaitbar=waitbar(0);
end

% determine which pdf to use
if isequal(type, 'THRU')
    % one phase is interesting for thru
    phases = mod(t_s,param.samples_per_ui);
    if phases==0, phases = param.samples_per_ui; end
else
    phases=1:samp_UI;
end

mxV = zeros(size(phases));
% we already found the phase in the PSD process for MMSE
if strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE
    if isequal(type, 'THRU')
        pdf=get_pdf_from_sampled_signal(vs(:,phases), param.levels, delta_y); %#ok<AGROW>
    else
        pdf=get_pdf_from_sampled_signal(vs(:,ixphase), param.levels, delta_y);
    end
else
    for k=phases
        pdf_samples(k)=get_pdf_from_sampled_signal(vs(:,k), param.levels, delta_y); %#ok<AGROW>
        mxV(k)=sqrt(sum( pdf_samples(k).x.^2.*pdf_samples(k).y)); % standard deviation of PDF
        progress = k/length(phases);
        if OP.DISPLAY_WINDOW, waitbar(progress, hwaitbar, ['processing COM pdf ' chdata.base ] ); figure(hwaitbar); drawnow; end
    end
    [UNUSED_OUTPU, pxi]=max(mxV); %#ok<ASGLU>
    pdf=pdf_samples(pxi);
end



if OP.DISPLAY_WINDOW
    close(hwaitbar);
end