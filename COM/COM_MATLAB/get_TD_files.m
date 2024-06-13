function [chdata, param] = get_TD_files(param, OP, num_fext, num_next, file_list)
% filename parsing and acquisition
%------------------------------------------------------------------
%----------put files names into chdata structure ---------
% The thru file has the index of 1
% crosstalk file are indexed from 2
% nxi is incremented each time a file is read in  so that nxi will end
filepath=[]; % path name for file
nxi=0; % file index
% get the THRU file
if size(file_list,2) ~= 0
    file_list(1)=strrep(file_list(1),'\', filesep); % OS file convention conversion
    [filepath, basename, fileext]=fileparts(file_list{1});
    
else
    if OP.RX_CALIBRATION == 1 || OP.FIXTURE_CALIBRATION
        h = msgbox('enter test channel s-parameter file'); set(h,'Color',[1 .85 0]);
        movegui(h,'northeast')
    end
        dir=fullfile(filepath, '*.csv');
        [basename,filepath]=uigetfile(dir,'input thru channel response .cav ');
        if filepath == 0
            error('No Thru file')
        end
    [UNUSED_OUTPUT, basename, fileext]=fileparts(basename); %#ok<ASGLU>
end
nxi=nxi+1;
chdata(nxi).filename = fullfile(filepath,  [basename fileext]);
chdata(nxi).ext = fileext;
[UNUSED_OUTPUT, dirname]=fileparts(filepath); %#ok<ASGLU>
%    chdata(nxi).base=[pth(max(strfind(pth,filesep))+1:end) '--' basename]; % add 1 directory back to basename
chdata(nxi).base=[OP.RUNTAG ' ' dirname  '--' basename]; % add 1 directory back to basename
% chdata(nxi).A=param.A_thru; % pam encoding amplitude reduction is do in reporting but is comprehending in crosstalk PDF
chdata(nxi).type='THRU';
chdata(nxi).ftr=param.fb*param.f_v;
param.base=chdata(nxi).base; %for print out function that don't pass chdata

% now get FEXT file names into chdata structure
kxi=nxi;
for nxi=kxi+1:num_fext+kxi
    lastfilepath=filepath;
    if size(file_list,2) ~= 0
        [filepath, basename, fileext]=fileparts(file_list{nxi});
    else
        if OP.RX_CALIBRATION == 1 || OP.FIXTURE_CALIBRATION
            h = msgbox('enter noise channel s-parameter file'); set(h,'Color',[1 .85 0]);
            movegui(h,'northeast')
        end
        if  param.tfx(1) == -1 && OP.ERL_ONLY && OP.ERL == 2;
            dir=fullfile(filepath, '*.csv');
            [basename,filepath]=uigetfile(dir,'input noise channel response .csv');
            if filepath==0
                error('Not enough NEXT files')
            end
        else
            dir=fullfile(filepath, '*.csv');
            if OP.RX_CALIBRATION == 1 || OP.FIXTURE_CALIBRATION
                [basename,filepath]=uigetfile(dir,'input noise channel response .csv');
            else
                [basename,filepath]=uigetfile(dir,['input fext channel response .csv #', num2str(nxi-kxi)]);
            end
            if filepath==0
                error('Not enough NEXT files')
            end
        end
        [UNUSED_OUTPUT, basename, fileext]=fileparts(basename); %#ok<ASGLU>
    end
    if isempty( filepath), filepath=lastfilepath; end
    chdata(nxi).filename = fullfile(filepath,  [basename fileext]);
    chdata(nxi).ext = fileext;
    [UNUSED_OUTPUT, dirname]=fileparts(filepath); %#ok<ASGLU>
    chdata(nxi).base=[OP.RUNTAG ' ' dirname '--' basename ];
    %     chdata(nxi).A=param.a_fext;
    chdata(nxi).ftr=param.fb*param.f_f;
    chdata(nxi).type='FEXT';
end
% now get NEXT file names into chdata structure
kxi=num_fext+kxi;
for nxi=kxi+1:num_next+kxi
    lastfilepath=filepath;
    if size(file_list,2) ~= 0
        [filepath, basename, fileext]=fileparts(file_list{nxi});
    else
        dir=fullfile(filepath, '*.csv');
        [basename,filepath]=uigetfile(dir,['input next channel response .csv ', num2str(nxi-kxi)]);
        if filepath==0
            error('Not enough NEXT files')
        end
        [UNUSED_OUTPUT, basename, fileext]=fileparts(basename); %#ok<ASGLU>
    end
    if isempty( filepath), filepath=lastfilepath; end
    chdata(nxi).filename = fullfile(filepath,  [basename fileext]);
    chdata(nxi).ext = fileext;
    [UNUSED_OUTPUT, dirname]=fileparts(filepath); %#ok<ASGLU>
    chdata(nxi).base=[OP.RUNTAG ' ' dirname '--' basename ];
    %     chdata(nxi).A=param.A_next;
    chdata(nxi).ftr=param.fb*param.f_n;
    chdata(nxi).type='NEXT';
end