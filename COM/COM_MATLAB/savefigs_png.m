function [ h ] = savefigs_png( param, OP )

%% find the figures
hw = waitbar(0,'Saving figures...');
h = findobj(0, 'Type', 'figure');
for ii=1:length(h)
    
    figname= get(h(ii), 'Name'); % use the figure name as file name
    if isempty(strfind(figname,param.base))
        figname = [figname ' ' OP.RUNTAG ' ' param.base ]; %#ok<AGROW>
    end
    if verLessThan('matlab', '8.4.0')
        figname = ['f_' num2str(h(ii)) '_' figname]; %#ok<AGROW>
    else
        figname = ['f_' num2str(h(ii).Number) '_' figname]; %#ok<AGROW>
    end
    figname = strrep(figname,':','-');
    figname = strrep(figname,' ','_');
    if OP.SAVE_FIGURES==1
        saveas(h(ii), fullfile(OP.RESULT_DIR, [figname '.png']));
    end
    %% get x y data
    if OP.SAVE_FIGURE_to_CSV==1
        h_L = findobj(h(ii),'Type','line'); % find handles to all the lines
        M=[]; %ncol=1;
        for nk=1:length(h_L)
            % get x and data for a line.
            x_data=get(h_L(nk),'xdata')';
            y_data=get(h_L(nk),'ydata')';
            % .........>> need to get data in the line structure (legend or label) for headers
            M=[M; x_data; y_data]; %#ok<AGROW>
        end
        csvwrite([OP.RESULT_DIR figname '.csv'],M);
        %      clear M y x header h_L
    end
    waitbar(ii/length(h),hw)
    
end

close(hw)
