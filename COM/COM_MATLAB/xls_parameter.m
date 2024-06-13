function p=xls_parameter(param_sheet, param_name, eval_if_string, default_value)
% helper function to read parameter values from XLS file. Uses names to find values.
if nargin<3, eval_if_string=0; end
[row, col]=find(strcmpi(param_sheet, param_name)); % RIM 08-26-2020 make case insensitive
if numel(row)*numel(col)==0
    if nargin<4
        missingParameter(param_name);
    else
        p = default_value;
    end
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
if ischar(p) && eval_if_string
    p=eval(p);
end
OP.SAVE_KEYWORD_FILE=0; % OP not passed ... maybe later set to 1 manually to get keyword file.
if OP.SAVE_KEYWORD_FILE

    if nargin<3 ||  ~exist('default_value','var') 
        default_value=p;
    end
    if isempty(default_value)
        default_value='-';
    end
    %%
    % Get call-stack info:
    stDebug = dbstack;
    callerFileName = stDebug(2).file;
    callerLineNumber = stDebug(2).line;
    % Open caller file:
    fCaller = fopen(callerFileName);
    % Iterate through lines to get to desired line number:
    for iLine = 1 : callerLineNumber
        % Read current line of text:
        currLine = fgetl(fCaller);
    end
    % (currLine) now reflects calling desired code: display this code:
%     fprintf('Complete text of calling code is : ''%s''\n',currLine);
    % Close caller file:
    left_side=currLine(1:strfind(currLine,'=')-1);
    cmt_side=currLine(strfind(currLine,'%')+1:end);
    if isempty(cmt_side), cmt_side=' ';end
    fclose(fCaller);
    
    if ~ischar(default_value)
        default_str=sprintf('%g ',default_value);
    else
        default_str=default_value;
    end
    if ~isfile('keyworklog.mat')
        save_p=param_name;
        save_d=default_str;
        save_r=left_side;
        save_c=cmt_side;
        param_name = {'keyword'};
        default_str =  {'default'};
        left_side={'matlab variable'};
        cmt_side={'info'};
        save('keyworklog.mat','left_side', 'param_name','default_str','cmt_side');
        param_name=save_p;
        default_str=save_d;
        left_side=save_r;
        cmt_side=save_c;
        data=load('keyworklog.mat');
    else
        load('keyworklog.mat');
    end
    data.left_side = [ data.left_side; left_side];
    data.param_name = [data.param_name; param_name];
    data.default_str = [data.default_str; default_str ];
    data.cmt_side = [ data.cmt_side; cmt_side];
    if length(data.default_str)~=length(data.default_str)
        a=1;
    end
    T=table(data.left_side, data.param_name, data.default_str, data.cmt_side);
    save('keyworklog.mat','data');
    writetable(T,[ 'keywords_' date '.csv' ]);
end