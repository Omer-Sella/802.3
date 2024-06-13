function Bathtub_Contribution_Wrapper(COM_SNR_Struct,Noise_Struct,param,chdata,OP)

% display bathtub curves in one axis per test case.
case_number=param.package_testcase_i;
if ~OP.COM_CONTRIBUTION_CURVES
    figure_name =  'Voltage bathtub curves';
    fig=findobj('Name', figure_name);
    if isempty(fig), fig=figure('Name', figure_name); end
    figure(fig);set(gcf,'Tag','COM');
    movegui(fig,'south')
    hax = subplot(length(OP.pkg_len_select), 1, case_number);
    plot_bathtub_curves( hax ...
        , COM_SNR_Struct.A_s ...
        , Noise_Struct.sci_pdf ...
        , Noise_Struct.cci_pdf ...
        , Noise_Struct.isi_and_xtalk_pdf ...
        , Noise_Struct.noise_pdf ...
        , Noise_Struct.jitt_pdf ...
        , COM_SNR_Struct.combined_interference_and_noise_pdf ...
        , param.delta_y ...
        );
    set(hax, 'tag', 'BTC');
    title(hax, sprintf('case %d VBC: %s ', case_number, regexprep([chdata(1).base,' '],'_',' ')));
    ylim(hax, [param.specBER/10 1]);
    % show BER target line
    hp=plot(get(hax, 'xlim'), param.specBER*[1 1], 'r:');
    set(get(get(hp,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
else
    figure_name =  'COM Contributions (Rough Allocations)';
    fig=findobj('Name', figure_name);
    if isempty(fig), fig=figure('Name', figure_name); end
    figure(fig);set(gcf,'Tag','COM');
    movegui(fig,'south')
    hax = subplot(length(OP.pkg_len_select), 1, case_number);
    
    plot_pie_com( hax ...
        , COM_SNR_Struct.A_s ...
        , Noise_Struct.sci_pdf ...
        , Noise_Struct.cci_pdf ...
        , Noise_Struct.isi_and_xtalk_pdf ...
        , Noise_Struct.noise_pdf ...
        , COM_SNR_Struct.combined_interference_and_noise_pdf ...
        , param.delta_y, param...
        );
    set(hax, 'tag', 'BTC');
    title(hax, sprintf('case %d rough COM impact: %s ', case_number, regexprep([chdata(1).base,' '],'_',' ')));
end

if OP.DEBUG && OP.DISPLAY_WINDOW && OP.RX_CALIBRATION==0
    btc_axes = findobj('tag', 'BTC');
    if ~isempty(btc_axes), linkaxes(btc_axes, 'x'); end
end