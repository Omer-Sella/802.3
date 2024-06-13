function half_UI=get_center_of_UI(samples_per_UI)

%half_UI reveals which value to use for the center of the UI.  For eye
%width calculations, it is necessary to place the cursor in the center of the
%UI window to ensure a 0 crossing on both left/right inside the window.
%This function was written in order to support even and odd samples_per_UI
%and to prevent the ambiguity of using samples_per_UI/2 vs. samples_per_UI/2+1

%The UI window goes from 0 to 1 with 1/samples_per_UI steps
UI_window=0:1/samples_per_UI:1-1/samples_per_UI;
%the center of the UI is sample closest to 0.5
[temp_diff,half_UI]=min(abs(UI_window-0.5));
function results= get_cm_noise(M,PR,L,BER,OP)

if ~exist('OP')
    OP.DC_norm_test=0;
    OP.DISPLAY_WINDOW=1;
end
param.BinSize=1e-5;
PR_test=-inf;
PR_fom_best=-inf;
% hwaitbar=waitbar(0);
for ki=1:M
    progress = ki/M;
    % if OP.DISPLAY_WINDOW
    %     waitbar(progress, hwaitbar, 'DM to CM computing'); figure(hwaitbar); drawnow;
    % else
    %     if ~mod(progress*100,1), fprintf('%i%% ', progress*100 );end
    % end
    tps=PR(ki:M:end);
    if OP.DC_norm_test
        PR_fom=(norm(tps));
    else
        testpdf=get_pdf_from_sampled_signal( tps,L, param.BinSize*10 );
        cdf_test=cumsum(testpdf.y);
        PRn_test=(-testpdf.x(find(cdf_test>=BER,1,'first')));
        PR_fom=PRn_test;
    end
    if PR_fom > PR_fom_best
        PR_fom_best=PR_fom;
        best_ki=ki;
    end
    if ~OP.DC_norm_test
        results.DCn=PR_fom_best;
        results.DCn_pdf=testpdf;
        results.DCn_cdf=cdf_test;
    else
        results.DCn=PR_fom_best;
    end
    results.DCn_p2p=max(PR)-min(PR);
end
