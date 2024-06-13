function [p,found]=xls_parameter_txffe(param_sheet, param_name)
% pretty much the same as "xls_parameter" but this is only used to dynamically find txffe
% to make the search dynamic, the "found" output is returned to let the calling function know to stop searching

found=1;
[row, col]=find(strcmpi(param_sheet, param_name)); % RIM 08-26-2020 make case insensitive
if numel(row)*numel(col)==0
    p = 0;
    found=0;
elseif numel(row)*numel(col)>1
    % if there are several occurrences, use the first, but warn
    %     warning('COM:XLS_parameter:MultipleOccurrence', ...
    %         '%d occurrences of "%s" found. Using the first', numel(row), param_name);
    error('COM:XLS_parameter:MultipleOccurrence', ...
        '%d occurrences of "%s" found. Please recheck spreadsheet', numel(row), param_name);% RIM 01-0 8-20
    p=param_sheet{row(1), col(1)+1};
else
    p=param_sheet{row, col+1};
end
if ischar(p)
    p=eval(p);
end
