function results=com_ieee8023_93a(varargin)
% This is NOT an official IEEE document.
%% Implementation example of annex 93A IEEE Std 802.3™
%   http://www.ieee802.org/3/ck/public/adhoc/index.html
% result=com_ieee8023_93a(config_file, num_fext, num_next [, <s4p files>])
% - config_file: xls, xls, mat file which contains configuration settings
% - num_fext: number of FEXT s4p files in the listfigure(300+package_testcase_i);
% - num_next: number of NEXT s4p files  in the list
% - <s4p_files>: (1+num_fext+num_nefxt) file names. If not supplied, program
%   will ask for each of the files interactively.opupu
%
% This program is intended for the development of standard specifications
% and reflects activity of IEEE P802.3bj, .3by, .3bm, .3bs, .3cd, .3ck
% found in Annex 93A IEEE Std 802.3™ and project =updates
% original proposal for COM may be found at
% http://www.ieee802.org/3/bj/public/jul12/mellitz_01_0712.pdf in July 2012
% from the following authors and affiliations in 2012.
%      Richard Mellitz, Intel Corporation
%      Charles Moore, Avago Technologies
%      Mike Dudek, QLogic Corporation
%      Mike Peng Li, Altera Corporation
%      Adee Ran, Intel Corporation
%
% Some of the authors and Contributors:
%   Adee Ran
%   Richard Mellitz
%   Yasuo Hidaka
%   John Ewen
%   Bill Kirkland
%   Adam Gregory
%   Howard Heck
%   Jingbo Li
%   Adam Healey
%   Matt Brown
%   Sameh Elnagar
%   Hossein Shakiba
zzz_list_of_changes()

%% Opening Preface
% acquire parsing command string and set up OP control structure. Then read in files
close(findall(0, 'tag', 'TMWWaitbar', '-or', 'tag', 'COM'));
try %  version number at end of call string
    cmdfile=mfilename;
    hindx=strfind(mfilename,'_');
    ver=cmdfile(hindx(end)+1:end);
    output_args.code_revision = [ver(1), '.',ver(2:end)];
catch
    output_args.code_revision ='';
end
teststr='';
OP.TESTING=0;
if OP.TESTING == 1 % set to 1 or pre release
    teststr='testing';
    testmsg=sprintf('Evaluation Copy: COM%s%s\n',output_args.code_revision,teststr);
    htest = msgbox(testmsg);
    set(htest,'Color','y', 'tag', 'COM');movegui(htest,'northeast');
end
display('This is NOT an official IEEE document.')
fprintf('Revision:<strong> %s%s </strong>This is a computation example for exploring COM and ERL  \n',output_args.code_revision,teststr)
disp(' for projects like IEEE P802.3bj/b/bs/cd/ck with some exploratory extensions and is not normative or official')
t0=tic;
set(0,'defaulttextinterpreter','none') % prevents subscripting character in displayed messages
% reset to tex on exit
%% file_setup
%%
% need to see what happens for version 8
if verLessThan('matlab', '7.4.1')
    error('Matlab version 7.4 or higher required')
end

results=[];

%% New Command Line parser
[config_file,num_fext,num_next,Remember_keyword,OP,varargin]=COM_CommandLine_Parse(OP,varargin{:});


%% get the first 3 arguments and allow for interactive input.
if isempty(config_file)
    config_file=input('Enter config XLS file or return will just pop a window to ask for the XLS file]:  ','s');
    if isempty(config_file)
        [config_file, config_file_path] = uigetfile([{ '*.xls;*.xlsx'} ; {'*.mat'}],'INPUT CONFIG FILE .xls');
    else
        [config_file_path,cname,cext]=fileparts(config_file);
        config_file=[cname cext];
    end
    if config_file==0
        % cancel - exit gracefully
        return;
    end
    config_file = fullfile(config_file_path, config_file);
end
output_args.config_file = config_file;
OP.SAVE_KEYWORD_FILE=0;
if OP.SAVE_KEYWORD_FILE
    if exist('keyworklog.mat','file')==2
        delete('keyworklog.mat');
    end
end
[param, OP] = read_ParamConfigFile(config_file,OP);
if OP.CONFIG2MAT_ONLY
    return;
end
if isempty(num_fext)
    if OP.RX_CALIBRATION
        num_fext=1;
        display('First prompt is for the measured test thru channel and following prompt is for Rx noise path channel')
    else
        if param.tfx(1) == -1 && OP.ERL_ONLY && OP.ERL == 2 % OP.ERL=2 is for package ERL, param.tfx = 1 means fixture time not set and need to be determinind for test fixture channel
            num_fext=1;
            display('First prompt is for the s2p measured data following prompt is for s4p of of the test fixtrure channel')
        elseif ~OP.ERL_ONLY
            num_fext=input('How many FEXT channels are to be entered? [return means no FEXT] ');
        else
            num_fext=0;
        end
    end
    if isempty(num_fext)==1, num_fext=0; end
end
if isempty(num_next)
    if param.tfx(1) == -1 && OP.ERL_ONLY && OP.ERL == 2 % OP.ERL=2 is for package ERL, param.tfx = 1 means fixture time not set and need to be determinind for test fixture channel
        num_next=0;
    else
        if OP.RX_CALIBRATION
            num_next=0;
        elseif ~OP.ERL_ONLY
            num_next=input('How many NEXT channels are to be entered? [return means no NEXT] ');
        else
            num_next=0;
        end
    end
    if isempty(num_next)==1, num_next=0; end
end
% Allow string inputs for running compiled version from OS command-line
if ischar(num_fext), num_fext=str2double(num_fext); end
if ischar(num_next), num_next=str2double(num_next); end
xtk=num_fext+num_next; % total number of crosstalk aggressors
param.num_next=num_next;
param.num_fext=num_fext;
param.num_s4p_files=num_fext+num_next+1;
% checking for data when running for rx compliance BBN calibration
if OP.RX_CALIBRATION == 1
    if num_fext ~=1
        h = msgbox('One and only noise path channel is required'); set(h,'Color',[1 .85 0]);
        movegui(h,'northwest')
        set(h,'Tag','COM') % RIM 06-13-2022 ... tak msg box for closing later
        if OP.DEBUG ~= 1
            return
        end
    end
    h = msgbox('Please make sure the measured "sigma_RJ", A_DD, and SNR_TX" fields in the config xls file have been modified from the Tx measurement. '); set(h,'Color',[0 1 1]);
    movegui(h,'southeast')
    set(h,'Tag','COM') % RIM 06-13-2022 ... tak msg box for closing later
end

if   param.tfx(1) == -1 && OP.ERL_ONLY && OP.ERL == 2
    if num_fext ~=1
        h = msgbox('One and only test channel is required'); set(h,'Color',[1 .85 0]);
         set(h,'Tag','COM') % RIM 06-13-2022 ... tak msg box for closing later
        movegui(h,'northwest')
        if OP.DEBUG ~= 1
            return
        end
    end
    h = msgbox('The test fixture file is use to gate measurements '); set(h,'Color',[0 1 1]);
    movegui(h,'southeast')
     set(h,'Tag','COM') % RIM 06-13-2022 ... tak msg box for closing later
end


% create result directory if needed
if ~exist(OP.RESULT_DIR,'dir'); mkdir(OP.RESULT_DIR); end
% allow finite impulse response input rather that s-parameters with
% OP.EXTERNAL = true. However the use_external_IR function is not provided
if ~isempty(varargin) % process case where file names are passed in function call
    if strfind(upper(char(varargin(1))),'EXTERNAL_IR') ~= 0
        error('External IR mode is no longer supported');
        %OP.EXTERNAL = true;
        %OP.GET_FD = 0;
        %ir1a= varargin(2);
        %ex_var = varargin(3);
        %[chdata OP param ]  = use_external_IR(param, OP ,num_fext,num_next,0,ir1a,ex_var);
    else
        if OP.TDMODE
            OP.GET_FD=false;
        end
        if length(varargin) < xtk +1 % check that number of varargin arguments passed is at least number of crosstalk files+1 (thru)
            error('files must include next + fext + a thru');
        end
        %% eveluate any extra arguments as possible modifications of parameters
        extra_args = varargin(xtk+2:end);
        for k=1:2:floor(length(extra_args)/2)*2
            try
                orig_value_is_str = 1;
                orig_value=eval(extra_args{k});
                if ~ischar(orig_value)
                    orig_value_is_str = 0;
                    orig_value=mat2str(orig_value);
                end
            catch eval_err
                if isequal(eval_err.identifier, 'MATLAB:nonExistentField')
                    % trying to modify a nonexistent parameter - probably a
                    % typo. save the user from his error.
                    error('COM:BadExtraParameter', 'Attempted override of a non-existing parameter %s.', extra_args{k});
                else
                    % unexpected condition
                    rethrow(eval_err);
                end
            end
            try
                if orig_value_is_str
                    mod_string = sprintf('%s = ''%s'';', extra_args{k}, extra_args{k+1});
                else
                    mod_string = sprintf('%s = %s;', extra_args{k}, extra_args{k+1});
                end
                eval(mod_string);
                fname=['mod_str' num2str(k)];
                % begin yasuo patch 2/11/2018
                % output_args.(fname)=mod_string;
                % If mod_string contains a comma, enclose it by double quotes to avoid misaligned column in the CSV output.
                
                % re-patch yasuo 3/18/2019
                % v2.56 if contains(mod_string,',')
                % v2.57 if isempty(strfind(mod_string,','))
                % Here, if-condition was inverted by the change of function from 'contains()' to 'isempty()'.
                % So, it is changed back by adding an '~' operator.
                % if isempty(strfind(mod_string,','))
                if ~isempty(strfind(mod_string,','))
                    output_args.(fname)=['"' mod_string '"'];
                else
                    output_args.(fname)=mod_string;
                end
                fprintf('Applied parameter modification: %s (override %s)\n', mod_string, orig_value);
            catch eval_err
                error(eval_err.identifier, 'Error evaluating "%s".', mod_string);
            end
        end
    end
end
%% Parameters computationally defined by values from the settings files
param.ui=1/param.fb;
param.sample_dt = param.ui/param.samples_per_ui;
param.sigma_X=sqrt( (param.levels^2-1)/ (3*(param.levels-1)^2) );
factor_3db=0.473037;
param.fb_BT_cutoff=factor_3db*param.f_r;
param.fb_BW_cutoff=param.f_r;
param.Tx_rd_sel=1;
param.Rx_rd_sel=2;
if isempty(param.snpPortsOrder) || any(isnan(param.snpPortsOrder))
    param.snpPortsOrder = [1 3 2 4]; % default order normally used.
end
%% size adjust vector parameters which may be entered as one element
param=parameter_size_adjustment(param,OP);

%% get input models
param.FLAG.S2P=0;
if OP.TDMODE
    OP.FIXTURE_CALIBRATION= 0;
    [chdata, param] = get_TD_files(param, OP, num_fext, num_next, varargin);
else
    OP.FIXTURE_CALIBRATION= 0;
    [chdata, param] = get_s4p_files(param, OP, num_fext, num_next, varargin);
    if any(strcmpi({chdata.ext},'.s2p'))
        param.FLAG.S2P=1;
    end
end

OP.SAVE_CMD_STR=1;
if OP.SAVE_CMD_STR
    cmd_str = save_cmd_line([Remember_keyword ''',''' config_file], chdata, num_fext,num_next,mfilename );
    setappdata(0,'cmd_str',cmd_str);
end
%% from here on, multiple package test cases are run. results will be saved separately.
results = cell(size(OP.pkg_len_select));
COM = inf;
min_COM=inf; % reset COM prior to calibration
% min_VEO = inf;
min_VEO_mV = inf;
max_VEC_dB = -inf;
threshold_DER=inf;
% begin yasuo patch 3/18/2019
threshold_DER_max = 0;	% reset worst DER
% end yasuo patch
sigma_bn=0;
DO_ONCE=true;
low_COM_found = 0;
% at this point only the impulse responses are needed. However vestiges of FD may be intermingled
while (OP.RX_CALIBRATION==1 || DO_ONCE==true)
    if ~DO_ONCE
        if abs(min_COM - param.pass_threshold)<0.1 || (sigma_bn==0 && min_COM < param.pass_threshold)
            break;
        elseif min_COM > param.pass_threshold
            % increase noise level linearly until low COM found; then perform binary search.
            if low_COM_found
                if OP.sigma_bn_STEP>0 % previous increase too small
                    OP.sigma_bn_STEP = OP.sigma_bn_STEP/2; % gearshift
                else % previously decrease too large
                    OP.sigma_bn_STEP = -OP.sigma_bn_STEP/2; % gearshift and change direction
                end
            end
        else % binary searchparam.Pkg_len_TX
            low_COM_found=1;
            if OP.sigma_bn_STEP>0 % previous increase too large
                OP.sigma_bn_STEP = -OP.sigma_bn_STEP/2; % gearshift and change direction
            else % previously decrease too small
                OP.sigma_bn_STEP = OP.sigma_bn_STEP/2; % gearshift
            end
        end
        min_COM = inf; % ignore previous iterations
        min_VEO_mV = inf;
        max_VEC_dB = -inf;
        sigma_bn = sigma_bn + OP.sigma_bn_STEP;
    end
    msgctr=1;
    for package_testcase_i = 1:length(OP.pkg_len_select)
        CSV_FILE=sprintf('%s%s_case%d_results.csv', OP.RESULT_DIR, chdata(1).base, package_testcase_i);
        package_testcase=OP.pkg_len_select(package_testcase_i);
        param.Pkg_len_TX = param.z_p_tx_cases(package_testcase,:);
        param.Pkg_len_NEXT = param.z_p_next_cases(package_testcase,:);
        param.Pkg_len_FEXT = param.z_p_fext_cases(package_testcase,:);
        param.Pkg_len_RX = param.z_p_rx_cases(package_testcase,:);
        param.AC_CM_RMS_TX= param.AC_CM_RMS(package_testcase);
        if param.PKG_Tx_FFE_preset ~=0
            param.Pkg_TXFFE_preset= param.PKG_Tx_FFE_preset(package_testcase,:);
        else
            param.Pkg_TXFFE_preset=0;
        end
        %         ki=package_testcase;
        % %         param.Pkg_Zc=[ param.pkg_Z_c(ki,1); param.pkg_Z_c(ki,2) ];
        %         param.Pkg_Zc=[ param.pkg_Z_c(ki,:) ];SDDp2p
        param.Pkg_Zc= param.pkg_Z_c;
        [cmele,centries] = size(param.Pkg_Zc);
        [mele, ncases] = size(param.Pkg_len_TX);
        if cmele ~=1 && centries ~=2 && mele ~= 1
            param.Pkg_Zc=reshape(param.Pkg_Zc,2,4);
        end
        param.package_testcase_i = package_testcase_i;
        
        %% Fill in chdata
        if OP.TDMODE
            [chdata, param ] = read_PR_files(param, OP, chdata);
            [chdata, param, SDDch, SDDp2p ] = TD_FD_fillin(param, OP, chdata); % fill in fd data to keep rest of SW happy
        else
            %fill in chada with s-parameters
            [chdata, SDDch, SDDp2p ] = read_s4p_files(param, OP, chdata);
            [chdata, param] = process_sxp(param, OP, chdata, SDDch);
        end
        if OP.BREAD_CRUMBS
            output_args.RL.f=chdata(1).faxis;  % RIM 07/19/2019 only use the first index
            output_args.RL.rl1=chdata(1).sdd11_raw; % RIM 07/19/2019 only use the first index
            if isfield(chdata(1),'sdd22_raw')% RIM 10/15/2019 only use the first index
                output_args.RL.r22=chdata(1).sdd22_raw; % RIM 07/19/2019 only use the first index
            end
            if isfield(chdata(1),'TX_RL')% RIM 10/09/2020 report Tx RL with RD
                output_args.RL.TXRL=chdata(1).TX_RL; %R IM 10/09/2020 report Tx RL with RD
            end
        end
        if param.FLAG.S2P, OP.ERL_ONLY =1;end
        
        %% Process TDR & ERL
        [output_args,ERL,min_ERL]=TDR_ERL_Processing(output_args,OP,package_testcase_i,chdata,param);
        if OP.ERL_ONLY
            results = cell(1);
            results{1} = output_args;
            rt=toc(t0);
            output_args.rtmin=rt/60;
            if 1
                fprintf('run time = %g min\n',output_args.rtmin)
            end
            if OP.CSV_REPORT ==1
                Write_CSV(output_args,CSV_FILE);
            end
            break;
        end
        
        %% FD processing s-parameter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % at this point sdd21 responses and faxis (frequency) array are defined
        %most operations now wrapped into FD_Processing function
        param.number_of_s4p_files=length(chdata);
        %ICN=0;
        output_args.ICN_mV=0;
        output_args.MDNEXT_ICN_92_46_mV=0;
        output_args.MDFEXT_ICN_92_47_mV=0;
        if OP.WC_PORTZ
            param.SNR_TX=param.SNDR(param.Tx_rd_sel);
        else
            param.SNR_TX=param.SNDR(package_testcase);
        end
        
        %TD Mode now also calls FD_Processing but skips the main parts
        [chdata,output_args]=FD_Processing(chdata,output_args,param,OP,SDDp2p,DO_ONCE);
        
        %% Convert from Frequency Domain to Time Domain
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if DO_ONCE
            if ~OP.TDMODE
                chdata=COM_FD_to_TD(chdata,param,OP);
                output_args.VMC_HF_mV=chdata(1).VCM_HF_struct.DCn*1000;
                output_args.SCMR_dB=chdata(1).SCMR;
            end
        end
        
        %% Determine equalization settings
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---------------------
        do_C2M=0;
        if param.T_O~=0 && param.Min_VEO_Test~=0
            do_C2M=1;
        end
        fom_result = optimize_fom(OP,param, chdata, sigma_bn,do_C2M);
        if fom_result.eq_failed ; return; end % RIM 12-20-2023
        
        %% Apply Equalization (returns pulse response with CTLE, TXLE, RXFFE)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_s=abs(fom_result.A_s); % this is the "s" in SNR (PAM4 gain in handled in last sections
        param.use_bmax=fom_result.best_bmax.';
        %AJG021820
        param.use_bmin=fom_result.best_bmin.';
        % Recommended Delta_y no larger than As/1000 or 0.01 mV
        param.current_ffegain=fom_result.best_current_ffegain;
        if OP.force_pdf_bin_size
            param.delta_y = OP.BinSize;
        else
            param.delta_y = min(A_s/1000, OP.BinSize);
        end
        % the pdf for PAM4 uses the full swing SBR but assigns voltage for the PDF accordingly
        if OP.RX_CALIBRATION, param.number_of_s4p_files=1; end
        
        chdata=Apply_EQ(param,fom_result,chdata,OP);
        
        PSD_results=[];
        if (strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
            OP.WO_TXFFE=1;
            %  chdata(1).eq_pulse_response includes rx FFE from Apply_EQ
            PSD_results = get_PSDs(PSD_results,chdata(1).eq_pulse_response,fom_result.t_s, fom_result.txffe,param.ctle_gdc_values(fom_result.ctle),param.g_DC_HP_values(fom_result.best_G_high_pass),param,chdata,OP);
            OP.WO_TXFFE=0;
            PSD_results = get_PSDs(PSD_results,chdata(1).eq_pulse_response,fom_result.t_s, fom_result.txffe,param.ctle_gdc_values(fom_result.ctle),param.g_DC_HP_values(fom_result.best_G_high_pass),param,chdata,OP);
        end
        %% Create ISI PDF & Individual Crosstalk PDFs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~OP.DISPLAY_WINDOW, fprintf('processing COM PDF '); end
        for i=1:param.number_of_s4p_files
            if ~OP.DISPLAY_WINDOW, fprintf('%d ', i); end
            if strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE

                pdf = get_pdf(chdata(i), param.delta_y, fom_result.t_s, param, OP,fom_result.PSD_results.iphase(i)) ;
            else
                pdf = get_pdf(chdata(i), param.delta_y, fom_result.t_s, param, OP,[]);
            end
            if OP.DEBUG && OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
                figure(150+package_testcase_i);set(gcf,'Tag','COM');
                subplot(2,1,2);
                pdf0.x=pdf.x(pdf.y~=0);
                pdf0.y=pdf.y(pdf.y~=0);
                semilogy(pdf0.x, pdf0.y,'disp', chdata(i).base);
                current_ylim=ylim; ylim([param.specBER/100, current_ylim(2)]);
                hold on; title('PDF')
                recolor_plots(gca);
            end
            
            chdata(i).pdfr=pdf;
            % reporting
            a=find(cumsum(chdata(i).pdfr.y) >1e-12,1,'first');
            chdata(i).maxquickpdf=(chdata(i).pdfr.y(a));
            
        end
        if ~OP.DISPLAY_WINDOW, fprintf('\n'); end
        
        %% Return final PDF & CDF and Package all noise parameters in Noise_Struct
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [PDF,CDF,Noise_Struct]=Create_Noise_PDF(A_s,param,fom_result,chdata,OP,sigma_bn,PSD_results);
        combined_interference_and_noise_pdf=PDF;
        combined_interference_and_noise_cdf=CDF;

              
        %% Calculate COM and other associated outputs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The noise and interference amplitude, A_ni, is the magnitude of the value of y0
        % that satisfies the relationship P(y0) = DER_0
        A_ni_ix=find(combined_interference_and_noise_cdf>param.specBER, 1, 'first');
        A_ni = abs(combined_interference_and_noise_pdf.x(A_ni_ix));

        % begin yasuo patch 3/18/2019
        % estimate DER at threshold COM
        threshold_ix=find(combined_interference_and_noise_pdf.x>-A_s/(10^(param.pass_threshold/20)),1);
        threshold_DER=combined_interference_and_noise_cdf(threshold_ix);
        threshold_DER_max = max(threshold_DER_max, threshold_DER);
        % end yasuo patch
        
        if OP.RX_CALIBRATION ==0 && OP.EW == 1
            [Left_EW,Right_EW,eye_contour,EH_T_C2M,EH_B_C2M]=COM_eye_width(chdata,param.delta_y,fom_result,param,OP,Noise_Struct,0);
            EW_UI=floor((Left_EW+Right_EW))/param.samples_for_C2M;
            if OP.DISPLAY_WINDOW && OP.DEBUG
                figure_name =  'Eye at DER0 estimate';
                fig=findobj('Name', figure_name);
                if isempty(fig), fig=figure('Name', figure_name); end
                figure(fig);set(gcf,'Tag','COM');
                movegui(fig,'southwest')
                plot(eye_contour)
                xlabel('UI %')
                ylabel('V')
            end
            
        else
            EW_UI=0;
            eye_contour=[];
        end
        if OP.MLSE==0
            if param.T_O ~=0
                eye_opening=EH_T_C2M-EH_B_C2M;
                A_ni=2*A_s-eye_opening;
                %eq 124E-4
                vec_arg=2*A_s/eye_opening;
                if vec_arg<eps
                    vec_arg=eps;
                end
                VEC_dB = 20*log10(vec_arg);
                COM=20*log10(2*A_s/A_ni);
                VEO_mV=eye_opening*1000;
                min_COM = min(min_COM, COM);
                min_VEO_mV = min(min_VEO_mV,VEO_mV);
                max_VEC_dB = max(max_VEC_dB, VEC_dB);
            else
                VEO_mV = 1000*(A_s-A_ni)*2;
                vec_arg=(A_s-A_ni)/A_s;
                if vec_arg<eps
                    vec_arg=eps;
                end
                VEC_dB = -20*log10(vec_arg);
                COM=20*log10(A_s/A_ni);
                min_COM = min(min_COM, COM);
                min_VEO_mV = min(min_VEO_mV,VEO_mV);
                max_VEC_dB = max(max_VEC_dB, VEC_dB);
            end
            MLSE_results=struct;
        else % MLSE case
            if OP.MLSE==1 
                [MLSE_results] =  MLSE(param,fom_result.DFE_taps(1),A_s,A_ni,PDF,CDF);
            elseif OP.MLSE==2
                [MLSE_results] =  MLSE_instu(param,fom_result.DFE_taps(1),A_s,A_ni,PDF,CDF);
            else
                warning('unsuported MLSE option')
            end
            if param.T_O ~=0
                eye_opening=EH_T_C2M-EH_B_C2M;
                A_ni=2*A_s-eye_opening;
                %eq 124E-4
                vec_arg=2*A_s/eye_opening;
                if vec_arg<eps
                    vec_arg=eps;
                end
                VEC_dB_orig = 20*log10(vec_arg); % was negative in 400 beta1 ... Fixed 2-2-23
                VEC_dB=MLSE_results.delta_com_CDF-20*log10( (10^(MLSE_results.delta_com_CDF/20)-1)*10^(VEC_dB_orig/20)+1)+VEC_dB_orig;
                COM_orig=20*log10(2*A_s/A_ni);
                COM=MLSE_results.COM_CDF;
                VEO_mV=eye_opening*1000;
                min_COM = min(min_COM, COM);
                min_VEO_mV = min(min_VEO_mV,VEO_mV);
                max_VEC_dB = max(max_VEC_dB, VEC_dB);
                output_args.delta_COM=MLSE_results.delta_com_CDF;
                output_args.delta_VEC=MLSE_results.delta_com_CDF-20*log10( (10^(MLSE_results.delta_com_CDF/20)-1)*10^(VEC_dB_orig/20)+1);
            else
                VEO_mV = 1000*(A_s-A_ni)*2;
                vec_arg=(A_s-A_ni)/A_s;
                if vec_arg<eps
                    vec_arg=eps;
                end
                VEC_dB_orig = -20*log10(vec_arg);
                VEC_dB=MLSE_results.delta_com_CDF-20*log10( (10^(MLSE_results.delta_com_CDF/20)-1)*10^(VEC_dB_orig/20)+1)+VEC_dB_orig;
                COM_orig=20*log10(A_s/A_ni);
                COM=MLSE_results.COM_CDF;
                min_COM = min(min_COM, COM);
                min_VEO_mV = min(min_VEO_mV,VEO_mV);
                max_VEC_dB = max(max_VEC_dB, VEC_dB);
                output_args.delta_COM=MLSE_results.delta_com_CDF;
            end
        end
        
        %% Create COM_SNR_Struct to hold the main COM outputs
        COM_SNR_Struct.A_s=A_s;
        COM_SNR_Struct.A_ni=A_ni;
        COM_SNR_Struct.threshold_DER=threshold_DER;
        COM_SNR_Struct.EW_UI=EW_UI;
        COM_SNR_Struct.COM=COM;
        COM_SNR_Struct.VEC_dB=VEC_dB;
        if OP.MLSE == 0
            COM_SNR_Struct.COM_orig=[];
            COM_SNR_Struct.VEC_dB_orig=[];
        else
            COM_SNR_Struct.COM_orig=COM_orig;
            COM_SNR_Struct.VEC_dB_orig=VEC_dB_orig;
        end
        COM_SNR_Struct.VEO_mV=VEO_mV;
        COM_SNR_Struct.combined_interference_and_noise_pdf=combined_interference_and_noise_pdf;
        COM_SNR_Struct.combined_interference_and_noise_cdf=combined_interference_and_noise_cdf;
        COM_SNR_Struct.eye_contour=eye_contour;
        
        
        %% Save TD
        if OP.SAVE_TD
            sbr=timeseries(fom_result.sbr,fom_result.t);
            if ~OP.TDMODE
                fir=timeseries(fom_result.IR,fom_result.t);
            end
            for i=1:param.number_of_s4p_files
                Pulses(i).uneq_responce= timeseries(chdata(i).uneq_pulse_response, chdata(i).t );
                Pulses(i).eq_responce= timeseries(chdata(i).eq_pulse_response, chdata(i).t );
                if ~OP.TDMODE
                    FIR(i).uneq_imp_response=  timeseries(chdata(i).uneq_imp_response, chdata(i).t );
                    FIR(i).eq_imp_response=  timeseries(chdata(i).eq_imp_response, chdata(i).t );
                end
            end
            if OP.TDMODE
                save( [OP.RESULT_DIR 'sbr_fir_' param.base '.mat'],'sbr','Pulses');
            else
                save( [OP.RESULT_DIR 'sbr_fir_' param.base '.mat'],'sbr','fir','Pulses','FIR')
            end
        end
        
        %% Bathtub/Contribution Plot
        if OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
            Bathtub_Contribution_Wrapper(COM_SNR_Struct,Noise_Struct,param,chdata,OP);
        end
        
        %% Msg management
        if ~exist('msg','var')
            msg=[];
        end
        if OP.DEBUG
            [ncases, mele]=size(param.z_p_tx_cases);
            switch param.flex
                case 4
                    msg = sprintf('%s: Case %g: z_p=(%g:%g:%g:%g, %g:%g:%g:%g, %g:%g:%g:%g, %g:%g:%g:%g) (TX, RX, NEXT, FEXT):\n', ...
                        msg,package_testcase_i, param.Pkg_len_TX, param.Pkg_len_RX, param.Pkg_len_NEXT, param.Pkg_len_FEXT ...
                        );
                case 2
                    msg = sprintf('%s: Case %g: z_p=(%g:%g, %g:%g, %g:%g, %g:%g) (TX, RX, NEXT, FEXT):\n', ...
                        msg,package_testcase_i, param.Pkg_len_TX(1:2), param.Pkg_len_RX(1:2), param.Pkg_len_NEXT(1:2), param.Pkg_len_FEXT(1:2) ...
                        );
                otherwise
                    msg = sprintf('%s: Case %g: z_p=(%g, %g, %g, %g) (TX, RX, NEXT, FEXT):', ...
                        msg, package_testcase_i, param.Pkg_len_TX, param.Pkg_len_RX, param.Pkg_len_NEXT, param.Pkg_len_FEXT ...
                        );
                    
            end
        else
            msg = sprintf('Case %d:', package_testcase_i );
        end
        
        if OP.TDMODE
            min_ERL=inf;
            ERL=[inf inf];
        end
        [msg] = end_display_control(msg,param,OP,output_args,COM,min_ERL,ERL, VEO_mV,VEC_dB,threshold_DER,OP.DISPLAY_WINDOW); % {} forces no ERL print
        
        
        %% Output Args
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        output_args=Output_Arg_Fill(output_args,sigma_bn,Noise_Struct,COM_SNR_Struct,param,chdata,fom_result,OP);
        rt=toc(t0);
        output_args.rtmin=rt/60;
        
        if OP.BREAD_CRUMBS
            output_args.OP=OP;
            output_args.param=param;
            output_args.chdata=chdata;
            output_args.fom_result = fom_result;
            output_args.PDF=PDF; % for exploration
            output_args.CDF=CDF; % for exploration
            output_args.MLSE_results=MLSE_results;
        end
        % results{package_testcase_i} = output_args;% moved RIM 04-14-2023
        
        %% making csv file
        if OP.CSV_REPORT ==1
            Write_CSV(output_args,CSV_FILE);
        end
        %% making mat file
        if(OP.DEBUG)
            save (sprintf('%s%s_case%d_results.mat', OP.RESULT_DIR, chdata(1).base, package_testcase_i), ...
                'output_args','param','OP');
        end
        if 1
            fprintf(' Die to die loss = %g dB \n',output_args.IL_db_die_to_die_at_Fnq)
            fprintf('run time = %g min \n',output_args.rtmin)
        end
        
        if nargout==0
            fprintf('<strong>--- Testcase %d results ---</strong>\n', package_testcase_i);
            disp(output_args)
        end
        
        if OP.BREAD_CRUMBS
            [my_path,rootname]=fileparts(chdata(1).filename);
            if ~isempty(OP.BREAD_CRUMBS_FIELDS)
                %Attempt to reduce the size of output_args.chdata by removing certain fields
                try
                    output_args.chdata=Bread_Crumb_Chdata_Reduction(output_args.chdata,OP.BREAD_CRUMBS_FIELDS);
                catch
                    fprintf('Failed to reduce output_args.chdata\n');
                end
            end
            save (sprintf('%s%s_case%d_results.mat', OP.RESULT_DIR, rootname, package_testcase), ...
                'output_args','param','OP');
        end
        
        results{package_testcase_i} = output_args; % moved to after chdata field reduction RIM 04-14-2023
    end
    [tmp] = end_display_control('WC All cases',param,OP,output_args,min_COM,min_ERL,ERL,min_VEO_mV,max_VEC_dB,threshold_DER,0);
    %%

    if OP.RX_CALIBRATION ==1
        sigma_hp= Noise_Struct.sigma_hp; % added for clause 162 else sigma_bn = sigma_hp (RIM 09-30-2022)
        display ([' LOOP with [sigma_bn sigma_hp] = [' num2str(sigma_bn) ' ' num2str(sigma_hp) '] performed with COM = ' num2str(min_COM) ])
    end
    DO_ONCE=false;
end

%% Final cleanup
if OP.DISPLAY_WINDOW
    savefigs(param, OP);
    set(0,'defaulttextinterpreter','tex'); % reset defaut text interpreter to tex
end

if OP.RX_CALIBRATION==1 % updated for clause 162 else sigma_bn = sigma_hp (RIM 09-30-2022)
    if ~param.f_hp==0
        fprintf ('Set Tx calibration noise(sigma_hp) rms voltage to %g mV\n', sigma_hp*1000);
        if OP.DISPLAY_WINDOW
            message=sprintf('Set Tx calibration noise (sigma_hp) rms voltage to %g mV.',sigma_hp*1000);
            hlast = msgbox(message,'sigma_hp','help');
            set(hlast,'Color','y', 'tag', 'COM');
        end
    else
        fprintf ('Set calibration noise (sigma_bn)rms voltage to %g mV\n', sigma_bn*1000);
        if OP.DISPLAY_WINDOW
            message=sprintf('Set calibration noise rms (sigma_bn) voltage to %g mV.',sigma_bn*1000);
            hlast = msgbox(message,'sigma_bn','help');
            set(hlast,'Color','y', 'tag', 'COM');
        end
    end
end

if length(results)==1, results = results{1}; end
redo_cmd_str=' redo string is: eval([''My_var_0 = '' getappdata(0,''cmd_str'')])';
disp(redo_cmd_str);
if isdeployed
    if OP.exit_if_deployed
        quit
    end
end
%%
%--------------------------------------------------------------------------
%--------------- Helper functions -----------------------------------------
%--------------------------------------------------------------------------
function chdata=Apply_EQ(param,fom_result,chdata,OP)

FB=param.fb;
FZ=param.CTLE_fz(fom_result.ctle);
FP1=param.CTLE_fp1(fom_result.ctle);
FP2=param.CTLE_fp2(fom_result.ctle);
GDC=param.ctle_gdc_values(fom_result.ctle);
if ~isempty(param.f_HP)
    FHP=param.f_HP(fom_result.best_G_high_pass);
end
if ~isempty(param.g_DC_HP_values)
    GDCHP=param.g_DC_HP_values(fom_result.best_G_high_pass);
end
if ~isempty(param.f_HP_Z)
    FHPZ=param.f_HP_Z(fom_result.ctle);
end
if ~isempty(param.f_HP_P)
    FHPP=param.f_HP_P(fom_result.ctle);
end
%Handle the scenario where the pulse response is not long enough to
%contain all DFE taps.  the SBR recorded in fom_result has the proper
%length
SBR_Len=length(fom_result.sbr);
if length(chdata(1).uneq_imp_response)<SBR_Len
    samples_added=SBR_Len-length(chdata(1).uneq_imp_response);
    chdata(1).uneq_imp_response(end+1:SBR_Len)=0;
    chdata(1).uneq_pulse_response(end+1:SBR_Len)=0;
    chdata(1).t(end+1:SBR_Len)=(1:samples_added)/param.fb/param.samples_per_ui+chdata(1).t(end);
end
for i=1:param.number_of_s4p_files
    % get quick PDF  results but only for THRU when in Rx calibration
    uneq_ir=chdata(i).uneq_imp_response;% includes packages, Hx, and Hr
    if OP.INCLUDE_CTLE==1
        switch param.CTLE_type
            case 'CL93'
                eq_ir = TD_CTLE(uneq_ir, FB, FZ, FP1, FP2, GDC, param.samples_per_ui);
            case 'CL120d'
                eq_ir = TD_CTLE(uneq_ir, FB, FZ, FP1, FP2, GDC, param.samples_per_ui);
                eq_ir = TD_CTLE(eq_ir, FB, FHP, FHP, 100e100 , GDCHP, param.samples_per_ui);
            case 'CL120e' % z has been adjusted  for gain
                eq_ir = TD_CTLE(uneq_ir, FB, FZ, FP1, FP2, GDC, param.samples_per_ui);
                eq_ir = TD_CTLE(eq_ir,FB, FHPZ, FHPP,1e99, 0, param.samples_per_ui);
        end
    else
        eq_ir=uneq_ir;
    end
    chdata(i).eq_imp_response=eq_ir;
    eq_pulse=filter(ones(1, param.samples_per_ui), 1, chdata(i).eq_imp_response);

    if isequal(chdata(i).type, 'FEXT') || isequal(chdata(i).type, 'THRU')
        eq_pulse = FFE( fom_result.txffe ,fom_result.cur-1 , param.samples_per_ui,  eq_pulse );
    end
    % chdata(i).ctle_imp_response
    if OP.RxFFE
        if isequal(upper(OP.FFE_OPT_METHOD),'MMSE')
            chdata(i).ctle_imp_response = FFE( fom_result.RxFFE ,fom_result.cur-1 , param.samples_per_ui,  eq_ir );
        end
        [ eq_pulse, C]=force(eq_pulse,param,OP,fom_result.t_s,fom_result.RxFFE);
    end
    chdata(i).eq_pulse_response=eq_pulse;% includes packages, Hf, and Hr, Ht, and Hffe(from tx)
end
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
function H_bt=Bessel_Thomson_Filter(param,f,use_BT)

if use_BT
    a = bessel( param.BTorder );
    acoef=fliplr( a );
    H_bt =a(1)./ polyval(acoef, (1i*f./(param.fb_BT_cutoff*param.fb)));
else
    H_bt=ones(1,length(f));
end
    
    
function chdata=Bread_Crumb_Chdata_Reduction(chdata,fields_file)

%This optional function reduces the size of output_args.chdata by parsing user supplied fields in a txt file
%The first line of the file must be #reduce or #include
%All subsequent lines are field names in chdata
%If using #reduce, the list of fields are the fields to remove from chdata
%If using #include, the list of fields are the fields to include in chdata
%
%Example file to remove the fields "sdd12_raw" "sdd21_raw" "sdd22_raw" "sdd11_raw"
%#reduce
%sdd12_raw
%sdd21_raw
%sdd22_raw
%sdd11_raw
%

fid=fopen(fields_file,'r');
file_data=textscan(fid,'%s','Delimiter','\n');

file_data=file_data{1};
fclose(fid);

%remove blank lines
L=cellfun('length',file_data);
file_data=file_data(L~=0);

%first line must be '#reduce' or '#include'
type=file_data{1};
field_names=file_data(2:end);
switch lower(type)
    case '#reduce'
        remove_fields=field_names;
    case '#include'
        all_fields=fieldnames(chdata);
        remove_fields=setdiff(all_fields,field_names);
    otherwise
        error('Bad first line.  Must be "#reduce" or "#include"');
end

%remove the "remove_fields" from chdata
for j=1:length(remove_fields)
    this_field=remove_fields{j};
    if isfield(chdata,this_field)
        chdata=rmfield(chdata,this_field);
    end
end
function [p_burst,p_error_propagation]=Burst_Probability_Calc(COM_SNR_Struct,DFE_taps,param,OP)

% an error burst of length N will cause each of the first N taps tap to mis-correct and create a PAM (2 or
% 4) noise term - depending on the N'th previous symbol, with double the tap voltage. From this we calculate
% the probability of staying in error state, i.e. burst of length N+1.

A_s=COM_SNR_Struct.A_s;
% initialize loop with uncorrelated noise and BER
error_propagation_noise_pdf{1}=COM_SNR_Struct.combined_interference_and_noise_pdf;  % PDF for burst of length 1 is the uncorrelated PDF

% Assume an error will occur if the noise excceds the available signal
% reduced by some dB. reduction is by COM threshold minus Error
% propagation margin (a positive EP margin reduces uncorrelated error probability
% below target BER).
error_threshold =  A_s./10^((param.pass_threshold-OP.COM_EP_margin)/20);
% Find the probability of this event by integration of the PDF. Use 1e-20 as a floor probabilty if noise PDF isn't wide enough.
x_error_propagation = find(error_propagation_noise_pdf{1}.x >= error_threshold, 1, 'first');
if isempty(x_error_propagation)
    p_error_propagation(1) = 1e-20;
else
    p_error_propagation(1) = sum(error_propagation_noise_pdf{1}.y(x_error_propagation:end));  % uncorrelated BER
end

sorted_abs_dfe_taps = sort(abs(DFE_taps), 'descend');
for k=2:min(param.ndfe, OP.nburst)
    % (arrays kept to allow tracking during development, though not really needed)
    if OP.use_simple_EP_model
        post_error_dfe_noise_pdf{k} = get_pdf_from_sampled_signal( 2*A_s*max(sorted_abs_dfe_taps), param.levels, param.delta_y ); %#ok<AGROW>
        error_propagation_noise_pdf{k} = conv_fct(error_propagation_noise_pdf{1}, post_error_dfe_noise_pdf{k}); %#ok<AGROW>
    else
        post_error_dfe_noise_pdf{k} = get_pdf_from_sampled_signal( 2*A_s*sorted_abs_dfe_taps(k-1), param.levels, param.delta_y ); %#ok<AGROW>
        error_propagation_noise_pdf{k} = conv_fct(error_propagation_noise_pdf{k-1}, post_error_dfe_noise_pdf{k}); %#ok<AGROW>
    end
    
    % Assume an error will propagate if this noise exceeds the threshold defined above
    x_error_propagation = find(error_propagation_noise_pdf{k}.x >= error_threshold, 1, 'first');
    if isempty(x_error_propagation)
        p_error_propagation(k) = 1e-20; %#ok<AGROW>
    else
        p_error_propagation(k) = sum(error_propagation_noise_pdf{k}.y(x_error_propagation:end)); %#ok<AGROW>
    end
end

% Assume an uncorrelated error will occur if the original noise exceeds
% the available signal reduced by pass_threhsold dB. Find the probability
% of this event by partial sum of the PDF.
% p_uncorrelated_error_i = find(combined_interference_and_noise_pdf.x >= A_s./10^(param.pass_threshold/20), 1, 'first');
% p_uncorrelated_error = sum(combined_interference_and_noise_pdf.y(p_uncorrelated_error_i:end));

% probability of bursts of different lengths
p_burst = cumprod(p_error_propagation);
function H_bw=Butterworth_Filter(param,f,use_BW)

if use_BW
    H_bw = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*f./(param.fb_BW_cutoff*param.fb));
else
    H_bw=ones(1,length(f));
end
function [CDF_ev] = CDF_ev(val,PDF,CDF)
index=find(PDF.x >= -val,1,'first');
CDF_ev=CDF(index);
function [CDF_inv_ev] = CDF_inv_ev(val,PDF,CDF)
index=find(CDF >= val,1,'first');
CDF_inv_ev=PDF.x(index);
function [config_file,num_fext,num_next,Remember_keyword,OP,varargin]=COM_CommandLine_Parse(OP,varargin)


keywords={'Legacy' 'TD' 'Config2Mat'};
Remember_keyword='Legacy';
OP.TDMODE=false;
OP.GET_FD=true;
OP.CONFIG2MAT_ONLY=false;
config_file='';
num_fext=[];
num_next=[];
if ~isempty(varargin)
    if ~ischar(varargin{1})
        error('First input must be a string');
    end
    keyword_idx=find(strcmpi(keywords,varargin{1}));
    if isempty(keyword_idx)
        %Legacy Mode
        [config_file,varargin]=varargin_extractor(varargin{:});
        [num_fext,varargin]=varargin_extractor(varargin{:});
        [num_next,varargin]=varargin_extractor(varargin{:});
    else
        %Keyword Mode
        my_keyword=varargin{1};
        Remember_keyword=my_keyword;
        varargin(1)=[];
        switch my_keyword
            
            case 'Legacy'
                [config_file,varargin]=varargin_extractor(varargin{:});
                [num_fext,varargin]=varargin_extractor(varargin{:});
                [num_next,varargin]=varargin_extractor(varargin{:});
            case 'TD'
                OP.TDMODE=true;
                OP.GET_FD=false;
                [config_file,varargin]=varargin_extractor(varargin{:});
                [num_fext,varargin]=varargin_extractor(varargin{:});
                [num_next,varargin]=varargin_extractor(varargin{:});
            case 'Config2Mat'
                OP.CONFIG2MAT_ONLY=true;
                [config_file,varargin]=varargin_extractor(varargin{:});
        end
    end
end
function chdata=COM_FD_to_TD(chdata,param,OP)

% get impulse responses which in interim step between equation for X(f) and
% H^(k)(t) without TX FFE or CTLE. These will we added later.
case_number=param.package_testcase_i;
for i=1:param.number_of_s4p_files
    %  RIM 2-01-2023 moved to FD_Processing
%     if OP.INCLUDE_FILTER % apply RX filtRaised_Cosine_Filterer
%         % Equation 93A-20 %%
%         %                     H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(i).faxis./(param.f_r*param.fb));
%         f=chdata(i).faxis;
%         %
%         H_bt=Bessel_Thomson_Filter(param,f,OP.Bessel_Thomson);
%         H_bw=Butterworth_Filter(param,f,OP.Butterworth);
%         H_RCos=Raised_Cosine_Filter(param,f,OP.Raised_Cosine); % conditionally include the RCos filter for all IR conversion using COM_FD_to_TD
%         H_txffe=      Tx_FFE_Filter(param,f,param.Pkg_TXFFE_preset); % RIM 08-18-2022 to add forced TX ffe per package case
%         H_r=H_bw.*H_bt.*H_RCos.*H_txffe; % RIM 08-18-2022 to add forced TX ffe per package case
%         chdata(i).sdd21=chdata(i).sdd21.*H_r;
%         if OP.DISPLAY_WINDOW
%             if i==1
%                 figure(300+param.package_testcase_i);
%                 subplot(3,1,1)
%                 plot(chdata(i).faxis/1e9, 20*log10(abs(squeeze(chdata(i).sdd21))), 'k-','linewidth',2, 'Disp','VTF (no Tx/Rx eq)')
%                 try
%                     legend('NumColumns',2)
%                     legend('location','south')
%                 catch
%                 end
%             end
%         end
%     end
    [chdata(i).uneq_imp_response, ...
        chdata(i).t, ...
        chdata(i).causality_correction_dB, ...
        chdata(i).truncation_dB] = s21_to_impulse_DC(chdata(i).sdd21 ,chdata(i).faxis, param.sample_dt, OP) ;
    if ~OP.RX_CALIBRATION || i==1 % DC (common to differentail model is not good used for RX_Calibrataion channel
        chdata(i).uneq_imp_response=chdata(i).uneq_imp_response*chdata(i).A; % adjust IRx for amplitude
        [chdata(i).uneq_CD_imp_response, ...
            chdata(i).t_DC, ...
            chdata(i).causality_correction_DC_dB, ...
            chdata(i).truncation__DC_dB] = s21_to_impulse_DC(chdata(i).sdc21 ,chdata(i).faxis, param.sample_dt, OP) ;
    end
    % adjust voltage derive here once it's decided what to use
    %------------------------------------------------------------
    % next find Pulse response (SBR) for each channel h^(k)(t)
    if ~OP.DISPLAY_WINDOW && i==1, fprintf('processing COM PDF '); end

    chdata(i).uneq_pulse_response=filter(ones(1, param.samples_per_ui), 1, chdata(i).uneq_imp_response);
    chdata(i).uneq_pulse_DC_response=filter(ones(1, param.samples_per_ui), 1, chdata(i).uneq_CD_imp_response);
    chdata(i).uneq_pulse_CD_response=chdata(i).uneq_pulse_DC_response*chdata(i).A;
    if 1 % not really CD but DC = DC if the channel already has a the Tx added
        % really need to add eq to the DC responce to calc rss.  This is a first pass estimate
        rss=-inf;
        for im=1:param.samples_per_ui
            rss=max(rss, norm( chdata(i).uneq_pulse_CD_response(im:param.samples_per_ui:end)));
        end
        chdata(i).CD_CM_RMS=rss*sqrt(param.sigma_X);
        chdata(i).VCM_HF_struct= get_cm_noise(param.samples_per_ui,chdata(i).uneq_pulse_CD_response,param.levels,param.specBER);
        chdata(i).SCMR=10*log10(max(chdata(1).uneq_pulse_response)^2/chdata(i).VCM_HF_struct.DCn^2);
    end
    if OP.DEBUG && OP.DISPLAY_WINDOW
        if OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
            figure(150+case_number);set(gcf,'Tag','COM');
            screen_size=get(0,'ScreenSize');
            pos = get(gcf, 'OuterPosition');
            set(gcf, 'OuterPosition', ...
                screen_size([3 4 3 4]).*[1 1 0 0] + pos([3 4 3 4]).*[-1 -1 1 1] ...
                - (case_number-1)*[0 20 0 0]);
            %movegui(gcf,'northeast')

            set(gcf, 'Name', sprintf('Case %d PR & PDF - %s', case_number, chdata(i).base));
            subplot(2,1,1);  hold on; % all plots on the same axes
            hp=plot(chdata(i).t, chdata(i).uneq_pulse_response,'Disp', chdata(i).base);
            hold on; % leave on for s-parameter problem finding. RIM 10-02-2023
            hp1=plot(chdata(i).t_DC, chdata(i).uneq_pulse_CD_response,'Disp', [ 'CD ' chdata(i).base ]) ;
        end
        % hide thru PR in order to show xtalk in a reasonable
        % scale. thru is shown in another plot.
        if isequal(chdata(i).type, 'THRU' ) && ~OP.RX_CALIBRATION %|| OP.RX_CALIBRATION % RIM 06-14-2022
            % set(hp, 'visible', 'off');
            % set(get(get(hp,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
        end
        title(sprintf('Unequalized Crosstalk and CD Conversion \n Pulse Responses'))
        ylabel('Volts')
        xlabel('seconds')

        recolor_plots(gca);
    else
        if param.ndfe~=0
            fprintf('%s\tUnequalized pulse peak = %.1f mV\n', chdata(i).base, 1000*max(abs(chdata(i).uneq_pulse_response)));
        end
    end

    fprintf('%s\tCausality correction = %.1f dB', chdata(i).base, chdata(i).causality_correction_dB);
    if OP.ENFORCE_CAUSALITY
        fprintf('\n');
    else
        fprintf(' (not applied)\n');
    end
    fprintf('%s\tTruncation ratio = %.1f dB\n', chdata(i).base, chdata(i).truncation_dB);

end
function [Left_EW,Right_EW,eye_contour,out_VT,out_VB]=COM_eye_width(chdata,delta_y,fom_result,param,OP,Struct_Noise,pdf_range_flag)


debug_plot=0;


samp_UI=param.samples_for_C2M;
half_UI=get_center_of_UI(samp_UI);
T_O=floor((param.T_O/1000)*param.samples_for_C2M);
start_sample=half_UI-T_O;
end_sample=half_UI+T_O;


%pdf_range is a placeholder to allow fractional UI calculation for optimize
%C2M speedup.  For regular COM_eye_width calls, pdf_range will be empty
%which will force pdf_range=1:samp_UI for stanard full UI calculation.
if pdf_range_flag
    pdf_range=[start_sample end_sample]; 
else
    pdf_range=[];
end

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




function [PDF,CDF,NS]=Create_Noise_PDF(A_s,param,fom_result,chdata,OP,sigma_bn,PSD_results)

%This block was originally in main COM function but was moved here for
%cleanup.  It returns the combined interference and noise PDF & CDF as well
%as a structure "NS" that contains all the noise parameters that are used
%in other places in COM

if OP.RX_CALIBRATION
    ctle_gain2 = (10.^(param.ctle_gdc_values(fom_result.ctle)/20) + 1i*chdata(2).faxis/param.CTLE_fz(fom_result.ctle)) ./ ...
        ((1+1i*chdata(2).faxis/param.CTLE_fp1(fom_result.ctle)).*(1+1i*chdata(2).faxis/param.CTLE_fp2(fom_result.ctle)));
    switch param.CTLE_type
        case 'CL93'
            H_low2=1;
        case 'CL120d' % this clause uses two gain indexes
            H_low2=(10.^(param.g_DC_HP_values(fom_result.best_G_high_pass)/20) +  1i*chdata(2).faxis/param.f_HP(fom_result.best_G_high_pass))./(1 + 1i*chdata(2).faxis/param.f_HP(fom_result.best_G_high_pass));
        case 'CL120e' % Z1 has been adjusted
            H_low2=(1 +  1i*chdata(2).faxis/f_HP_P(fom_result.ctle))./(1 + 1i*chdata(2).faxis/f_HP_Z(fom_result.ctle));
    end
    H_ctf2=H_low2.*ctle_gain2;
    [ sigma_ne, NS.sigma_hp] = get_sigma_noise( H_ctf2,  param, chdata, sigma_bn );
else
    sigma_ne=0;
end

NS.sigma_N = fom_result.sigma_N; % eta zero noise
if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
    if ~OP.SNR_TXwC0
        % Equation 93A-30 %% h0(ts0)= A_s*R_lm/(L-1)
        NS.sigma_TX = (param.levels-1)*A_s/param.R_LM*10^(-param.SNR_TX/20); % SNR_Tx and RLM
    else
        NS.sigma_TX = (param.levels-1)*A_s/fom_result.txffe(fom_result.cur) /param.R_LM*10^(-param.SNR_TX/20); % SNR_Tx from Adee
    end
else
    NS.sigma_TX =PSD_results.tn_rms;
end
% Equation 93A-41 %%
NS.sigma_G  = norm([param.sigma_RJ*param.sigma_X*norm(fom_result.h_J), NS.sigma_N, NS.sigma_TX]);
NS.sigma_rjit= param.sigma_RJ*param.sigma_X*norm(fom_result.h_J);

% Equation 93A-42 %%
% number of sigmas needed depends on the required BER.
if param.Noise_Crest_Factor == 0
    NS.ber_q = sqrt(2)*erfcinv(2*param.specBER);
else
    NS.ber_q=param.Noise_Crest_Factor;
end
NS.gaussian_noise_pdf = normal_dist(NS.sigma_G, NS.ber_q, param.delta_y);
% enable overriding the Q factor of the BBN instrument.
if OP.force_BBN_Q_factor
    NS.ne_noise_pdf = normal_dist(sigma_ne, OP.BBN_Q_factor, param.delta_y);
else
    NS.ne_noise_pdf = normal_dist(sigma_ne, NS.ber_q, param.delta_y);
end
NS.gaussian_noise_pdf = conv_fct(NS.gaussian_noise_pdf, NS.ne_noise_pdf);

% p_DD is computed using the procedure defined in 93A.1.7.1 with h(n)=A_DD*h_J(n)
NS.p_DD = get_pdf_from_sampled_signal(param.A_DD*fom_result.h_J, param.levels, param.delta_y);

% Equation 93A-43
NS.noise_pdf=conv_fct(NS.gaussian_noise_pdf, NS.p_DD);

gaussian_rjitt_pdf = normal_dist(NS.sigma_rjit, NS.ber_q, param.delta_y);
NS.jitt_pdf=conv_fct(gaussian_rjitt_pdf, NS.p_DD);

% Implementation of 93A.1.7.3 combination procedure
%  (effectively Equation 93A-44) %%

% Self-Channel Interference is thru residual result
NS.sci_pdf = chdata(1).pdfr;
sci_mxi=find(cumsum(NS.sci_pdf.y)>=param.specBER, 1, 'first');
NS.thru_peak_interference_at_BER=abs(NS.sci_pdf.x(sci_mxi));
sci_msi=find(cumsum(NS.sci_pdf.y)>=param.specBER, 1, 'first');
NS.sci_sigma=abs(NS.sci_pdf.x(sci_msi)/(erfcinv(2*param.specBER)*sqrt(2)));
if OP.RX_CALIBRATION ==0
    % Co-Channel Interference PDFs (for information only):
    % initialize to deltas
    MDNEXT_cci_pdf = d_cpdf(param.delta_y, 0, 1);
    MDFEXT_cci_pdf = d_cpdf(param.delta_y, 0, 1);
    % serially convolve FEXT/NEXT PDFs
    for k=2:param.number_of_s4p_files
        if isequal(chdata(k).type, 'NEXT')
            MDNEXT_cci_pdf = conv_fct(MDNEXT_cci_pdf, chdata(k).pdfr);
        else % ... must be FEXT
            MDFEXT_cci_pdf = conv_fct(MDFEXT_cci_pdf, chdata(k).pdfr);
        end
    end
    
    % find "peaks" of MDNEXT/MDFEXT for reporting
    mdnxi=find(cumsum(MDNEXT_cci_pdf.y)>=param.specBER, 1, 'first');
    NS.MDNEXT_peak_interference=abs(MDNEXT_cci_pdf.x(mdnxi));
    mdfxi=find(cumsum(MDFEXT_cci_pdf.y)>=param.specBER, 1, 'first');
    NS.MDFEXT_peak_interference=abs(MDFEXT_cci_pdf.x(mdfxi));
    
    % Combined crosstalk effect
    NS.cci_pdf = conv_fct(MDFEXT_cci_pdf, MDNEXT_cci_pdf);
    cci_mxi=find(cumsum(NS.cci_pdf.y)>=param.specBER, 1, 'first');
    cci_msi=find(cumsum(NS.cci_pdf.y)>=param.specBER, 1, 'first');
    NS.cci_sigma=abs(NS.cci_pdf.x(cci_msi)/(erfcinv(2*param.specBER)*sqrt(2)));
    NS.crosstalk_peak_interference_at_BER=abs(NS.cci_pdf.x(cci_mxi));
    % combine cci and sci
    NS.isi_and_xtalk_pdf = conv_fct(NS.sci_pdf, NS.cci_pdf);
else
    % for calibration there is no cci
    NS.isi_and_xtalk_pdf=NS.sci_pdf;
end

mxi=find(cumsum(NS.isi_and_xtalk_pdf.y)>=param.specBER, 1, 'first');
NS.peak_interference_at_BER=abs(NS.isi_and_xtalk_pdf.x(mxi));


% Equation 93A-45
combined_interference_and_noise_pdf = conv_fct(NS.isi_and_xtalk_pdf, NS.noise_pdf);
PDF=combined_interference_and_noise_pdf;

% Equation 93A-37
combined_interference_and_noise_cdf=cumsum(combined_interference_and_noise_pdf.y);
CDF=combined_interference_and_noise_cdf;
function [hctf] = FD_CTLE(freq, fb, f_z, f_p1, f_p2, kacdc_dB)
hctf = ( 10.^(kacdc_dB/20) + 1i*freq/f_z ) ./ ( (1+1i*freq/f_p1) .* (1+1i*freq/f_p2)) ;

function [chdata,output_args]=FD_Processing(chdata,output_args,param,OP,SDDp2p,DO_ONCE)
%This function calculates various frequency domain metrics
%Mainly IL_fit, FOM_ILD, ICN, ICN_Fext, and ICN_Next
db = @(x) 20*log10(abs(x));
package_testcase=OP.pkg_len_select(param.package_testcase_i);
if OP.WC_PORTZ
    A_thru = param.a_thru(param.Tx_rd_sel);
    A_fext = param.a_fext(param.Tx_rd_sel);
    A_next = param.a_next(param.Tx_rd_sel);
else
    A_thru = param.a_thru(package_testcase);
    A_fext = param.a_fext(package_testcase);
    A_next = param.a_next(package_testcase);
end
for i=1:param.number_of_s4p_files
    if isequal(chdata(i).type, 'THRU')
        chdata(i).A=A_thru;
        chdata(i).Aicn=A_thru;
    elseif isequal(chdata(i).type, 'FEXT')
        chdata(i).A=A_fext;
        chdata(i).Aicn=param.a_icn_fext;
    elseif isequal(chdata(i).type, 'NEXT')
        chdata(i).A=A_next;
        chdata(i).Aicn=param.a_icn_next;
    end
end
if OP.TDMODE
    for i=1:param.number_of_s4p_files           % freq delta for integration
        chdata(i).delta_f=chdata(i).faxis(11)-chdata(i).faxis(10);
    end
end
if ~DO_ONCE
    return;
end
%Any new output_args fields set in this function should be initialized here as empty
output_args.fitted_IL_dB_at_Fnq = [];
output_args.cable__assembley_loss=[];
output_args.loss_with_PCB=[];
output_args.VIP_to_VMP_IL_dB_at_Fnq=[];
output_args.IL_dB_channel_only_at_Fnq=[];
output_args.VTF_loss_dB_at_Fnq=[];
output_args.IL_db_die_to_die_at_Fnq=[];
output_args.FOM_TDILN=[];
output_args.TD_ILN=[];
output_args.FOM_RILN=[];
output_args.FOM_ILD=[];
%TD_Mode is just a pass through to set the empty values and return
if ~OP.GET_FD
    return;
end
case_number=param.package_testcase_i;
f2=param.f2;
f1=param.f1;
MDFEXT_ICN=0; MDNEXT_ICN=0;
for i=1:param.number_of_s4p_files
    if OP.INCLUDE_FILTER % apply RX filtRaised_Cosine_Filterer
        % Equation 93A-20 %%
        %                     H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(i).faxis./(param.f_r*param.fb));
        f=chdata(i).faxis;
        %
        H_bt=Bessel_Thomson_Filter(param,f,OP.Bessel_Thomson);
        H_bw=Butterworth_Filter(param,f,OP.Butterworth);
        H_RCos=Raised_Cosine_Filter(param,f,OP.Raised_Cosine); % conditionally include the RCos filter for all IR conversion using COM_FD_to_TD
        H_txffe=      Tx_FFE_Filter(param,f,param.Pkg_TXFFE_preset); % RIM 08-18-2022 to add forced TX ffe per package case
        H_r=H_bw.*H_bt.*H_RCos.*H_txffe; % RIM 08-18-2022 to add forced TX ffe per package case
        chdata(i).sdd21=chdata(i).sdd21.*H_r;
        if OP.DISPLAY_WINDOW
            if i==1
                figure(300+param.package_testcase_i);
                subplot(3,1,1)   
                hold on
                plot(chdata(i).faxis/1e9, 20*log10(abs(squeeze(chdata(i).sdd21))), 'k-','linewidth',2, 'Disp',sprintf('VTF (no Tx/Rx eq)')')
                try
                    legend('NumColumns',2)
                    legend('location','south')
                catch
                end
            end
        end
    end
end
for i=1:param.number_of_s4p_files
    if i == 2
        PSXT(1:length(chdata(i).sdd21f))=0;
        MDFEXT(1:length(chdata(i).sdd21f))=0;
        MDNEXT(1:length(chdata(i).sdd21f))=0;
    end
    a=find(chdata(i).faxis(:)>=f2,1,'first');% RIM 01-12-21
    if isempty(a)
        f2=chdata(i).faxis(end);
        index_f2=length(chdata(i).faxis);
    else
        index_f2=a(1);
    end
    b=find(chdata(i).faxis(:)<=f1,1,'last');% RIM 01-12-21
    if isempty(b)
        f1=chdata(i).faxis(1);
        index_f1=1;
    else
        index_f1=b(1);
    end
    % R is the frequency dependent parameter for the sinc function use in the
    % PWF for ICN
    temp_angle=(param.samples_per_ui*param.sample_dt)*pi.*chdata(i).faxis;
    if(chdata(i).faxis(1)==0)
        temp_angle(1)=1e-20;% we don't want to divide by zero
    end
    SINC = sin(temp_angle)./temp_angle;
    PWF_data=SINC.^2;
    PWF_trf=(1+(chdata(i).faxis/chdata(i).ftr).^4).^-1;
    %// bw1=2.613126; bw2=3.4142136; bw3=2.613126;
    fr=param.f_r*param.fb;
    PWF_rx=(1+(chdata(i).faxis/fr).^8).^-1;
    PWF_highpass=1;
    % Equation 93A-57 %
    PWF=PWF_data.*PWF_trf.*PWF_rx.*PWF_highpass; % power weight function
    % freq delta for integration
    chdata(i).delta_f=chdata(i).faxis(11)-chdata(i).faxis(10);
    % from ba spec, this is basically ICN
    faxis_GHz = chdata(i).faxis/1e9;
    if isequal(chdata(i).type, 'THRU')
                [ILD_magft chdata(i).fit_orig] = get_ILN(chdata(i).sdd21f(index_f1:index_f2), chdata(i).faxis(index_f1:index_f2));
        % find fitted loss values by interpolation using full data, no indexing - Adee 2022-08-28
        [~, chdata(i).fit_orig] = get_ILN(chdata(i).sdd21f, chdata(i).faxis);
        fit_loss = interp1(chdata(i).faxis, -chdata(i).fit_orig, 1/param.ui/2);
        chdata(i).fit_ILatNq = fit_loss;
        output_args.fitted_IL_dB_at_Fnq = fit_loss;
        IL_interp = interp1(chdata(i).faxis, -20*log10(abs(chdata(i).sdd21f)), 1/param.ui/2);
        chdata(i).ILatNq = IL_interp;
        if OP.include_pcb
            cable_loss = interp1(chdata(i).faxis, -20*log10(abs(chdata(i).sdd21_orig)), 1/param.ui/2);
            loss_with_PCB = interp1(chdata(i).faxis, -20*log10(abs(chdata(i).sdd21_raw)), 1/param.ui/2);
            output_args.cable__assembley_loss=cable_loss;
            output_args.loss_with_PCB=loss_with_PCB;
        end
        Nq_loss=chdata(i).ILatNq;
        output_args.IL_dB_channel_only_at_Fnq=Nq_loss;
        % time domain ref RR = complex fit pulse
        if OP.COMPUTE_TDILN || OP.COMPUTE_RILN
            [ILD chdata(i).fit, TD_ILN ] = get_ILN_cmp_td(chdata(i).sdd21f(index_f1:index_f2), chdata(i).faxis(index_f1:index_f2),OP,param,chdata(i).A);
            FOM_TDILN = TD_ILN.SNR_ISI_FOM_PDF;
            FOM_ILN_complex= TD_ILN.FOM;
        end
        if OP.COMPUTE_TDILN || OP.COMPUTE_RILN
            [ILD chdata(i).fit, TD_ILN ] = get_ILN_cmp_td(chdata(i).sdd21f(index_f1:index_f2), chdata(i).faxis(index_f1:index_f2),OP,param,chdata(i).A);
            FOM_TDILN = TD_ILN.SNR_ISI_FOM_PDF;
            FOM_ILN_complex= TD_ILN.FOM;
        end
        if OP.COMPUTE_TDILN
            output_args.FOM_TDILN=FOM_TDILN;
            output_args.TD_ILN=TD_ILN; % struct
        end
        if OP.COMPUTE_RILN
            % Get RIL, RILN, and TD_RILN
            [RIL_struct]= capture_RIL_RILN(chdata);
            FOM_RILN=sqrt(chdata(i).delta_f/(param.f2-param.f1)*sum( PWF(index_f1:index_f2-1).*RIL_struct.RILN_dB(index_f1:index_f2-1)'.^2));
            output_args.FOM_RILN=FOM_RILN;           
            %---start. plotting ILN based on ILD and RILN % Hansel 10/18/2021
            plot_tdomain_debug= 0; % must have OP.COMPUTE_TDILN = 1 to use
            if plot_tdomain_debug== 1
                figure(988); set(gcf,'Tag','COM')
                ax_1= subplot(3,1,1);
                plot(TD_ILN.REF.t*1e9, TD_ILN.REF.PR*1e3,'disp','ref');
                hold on;
                plot(TD_ILN.FIT.t*1e9, TD_ILN.FIT.PR*1e3,'disp','fit');
                hold on;
                plot(TD_ILN.t*1e9, TD_ILN.ILN*1e3, 'k','disp','ref - fit (IL noise)');
                ylim([min(TD_RILN.ILN) max(TD_RILN.ILN)]*1e3);
                grid on;
                box on;
                legend('REF', 'FIT', 'TD\_ILN: ref - fit (IL noise)');
                xlabel('Time [nsec]');
                ylabel('Pulse Response [mV]');
                
                ax_2= subplot(3,1,2);
                plot(TD_RILN.REF.t*1e9, TD_RILN.REF.PR*1e3);
                hold on;
                plot(TD_RILN.t*1e9, TD_RILN.ILN*1e3, 'r');
                ylim([min(TD_RILN.ILN) max(TD_RILN.ILN)]*1e3);
                grid on;
                box on;
                legend('REF', 'TD\_RILN');
                xlabel('Time [nsec]');
                ylabel('Pulse Response [mV]');
                ax_3= subplot(3,1,3);
                plot(TD_ILN.t*1e9, TD_ILN.ILN*1e3, 'k');
                hold on;
                plot(TD_RILN.t*1e9, TD_RILN.ILN*1e3, 'r');
                ylim([min(TD_RILN.ILN) max(TD_RILN.ILN)]*1e3);
                grid on;
                box on;
                legend( 'TD\_ILN: ref - fit (IL noise)', 'TD\_RILN');
                xlabel('Time [nsec]');
                ylabel('Pulse Response [mV]');
                
                linkaxes([ax_1, ax_2, ax_3], 'x');
                ax_1.XLim = [0 max(TD_RILN.t)*1e9 ];
            end
            %---end. plotting ILN based on ILD and RILN
        end
        % Equation 93A-56 %
        FOM_ILD=sqrt(chdata(i).delta_f/(param.f2-param.f1)*sum( PWF(index_f1:index_f2).*ILD_magft.^2));
        output_args.FOM_ILD=FOM_ILD;
        if OP.DEBUG
            if OP.DISPLAY_WINDOW                
                figure(300+case_number);
                set(gcf,'Tag','COM')
                screen_size=get(0,'ScreenSize');
                pos = get(gcf, 'OuterPosition');
                set(gcf, 'Name', [sprintf('%.3gdB IL Channel: ',Nq_loss) 'Raw frequency-domain data'], 'OuterPosition', ...
                    screen_size([3 4 3 4]).*[0 1 0 0] + pos([3 4 3 4]).*[0 -2 1 2] ...IL fit
                    - (case_number-1)*[0 20 0 0]);
                subplot(3,1,1)
                title('Losses')
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd21f))), 'b', 'LineWidth', 3, 'Disp','IL passed s-params')
                hold on
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd21p_nodie))), 'm-', 'Disp','IL die to die (w pkg/brds)')
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd21p))), 'b-', 'Disp','IL dB VIP to VMP')               
                ylim(get(gca, 'ylim'));
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd11))),'c','Disp','RL11')
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd22))),'m','Disp','RL22')
                subplot(3,1,3)
                plot(faxis_GHz(index_f1:index_f2), ILD_magft,'Disp','ILD')
                if OP.PLOT_CM
                    if case_number ==1
                        h350=figure(350);set(gcf,'Tag','COM')
                        screen_size=get(0,'ScreenSize');
                        pos = get(gcf, 'OuterPosition');
                        set(gcf, 'OuterPosition',screen_size([3 4 3 4]).*[1  1 0 0] + pos([3 4 3 4]).*[-2 -2 2 1])
                        movegui(gcf,'center');
                        htabgroup350 = uitabgroup(h350);
                        htab1 = uitab(htabgroup350, 'Title', 'CM Through Losses');
                        hax1 = axes('Parent', htab1);
                        set(h350,'CurrentAxes',hax1)
                        hold on
                        set(gcf,'Tag','COM')
                        screen_size=get(0,'ScreenSize');
                        pos = get(gcf, 'OuterPosition');
                        title('IL & CM Losses')
                        base=strrep(chdata(i).base,'_',' ');
                        plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd21f))), 'b', 'LineWidth', 3, 'Disp',  [ 'sdd21(IL) TP0-TP5 ' base])
                        plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdc21_raw))),   'LineWidth', 2, 'Disp',[ 'sdc21 TP0-TP5 '  base])
                        ylabel('dB')
                        xlabel('GHz')
                        legend show
                        legend('Location','eastoutside')
                        hold on
                        grid on
                        if param.number_of_s4p_files > 1
                            htab2 = uitab(htabgroup350, 'Title', 'CM Crosstalk Losses');
                            hax2 = axes('Parent', htab2);
                            htab3 = uitab(htabgroup350, 'Title', 'Crosstalk Losses');
                            hax3 = axes('Parent', htab3);
                        end
                        
                    end
                end
            else
                display(['Insertion Loss at Nyquist = ', num2str(chdata(i).ILatNq)])
            end
        end
    else % NEXT or FEXT
        if isequal(chdata(i).type, 'FEXT')
            MDFEXT=sqrt(abs(chdata(i).sdd21f).^2+MDFEXT.^2); % power sum xtk
            MDFEXT_ICN=sqrt(2*chdata(i).delta_f/param.f2*sum( chdata(i).Aicn^2*PWF(index_f1:index_f2).*abs(MDFEXT(index_f1:index_f2)).^2)); %eq 46
            output_args.MDFEXT_ICN_92_47_mV=MDFEXT_ICN*1000;
        elseif isequal(chdata(i).type, 'NEXT')
            MDNEXT=sqrt(abs(chdata(i).sdd21f).^2+MDNEXT.^2); % power sum xtk
            MDNEXT_ICN=sqrt(2*chdata(i).delta_f/param.f2*sum( chdata(i).Aicn^2*PWF(index_f1:index_f2).*abs(MDNEXT(index_f1:index_f2)).^2)); %eq 47
            output_args.MDNEXT_ICN_92_46_mV=MDNEXT_ICN*1000;
        end
        PSXT=sqrt((abs(chdata(i).sdd21f)*chdata(i).Aicn).^2+PSXT.^2); % power sum xtk
        ICN=sqrt(2*chdata(i).delta_f/param.f2*sum( PWF(index_f1:index_f2).*abs(PSXT(index_f1:index_f2)).^2));
        output_args.ICN_mV=ICN*1000;
        ICN_test=norm([MDFEXT_ICN MDNEXT_ICN]);
        if  OP.PLOT_CM && OP.DISPLAY_WINDOW
            if case_number ==1
                %                             htab2 = uitab(htabgroup350, 'Title', 'CM Crosstalk Losses');
                %                             hax2 = axes('Parent', htab2);
                set(h350,'CurrentAxes',hax2)
                hold on
                title('CM Losses')
                base=strrep(chdata(i).base,'_',' ');
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdc21_raw))), 'LineWidth', 2, 'Disp',[ 'sdc21 TP0-TP5 ' base ])
                legend('Location','eastoutside')
                hold on
                grid on
                set(h350,'CurrentAxes',hax3)
                plot(faxis_GHz, 20*log10(abs(squeeze(chdata(i).sdd21_raw))), 'LineWidth', 2, 'Disp',[ 'sdd21 TP0-TP5 ' base ])
                legend('Location','eastoutside')
                hold on
                grid on
            end
        end
    end
end  % for loop
ICN_test=norm([MDFEXT_ICN MDNEXT_ICN]);
if OP.DEBUG && OP.DISPLAY_WINDOW
    figure(300+case_number);set(gcf,'Tag','COM');
    if param.number_of_s4p_files > 1
        scale=1/chdata(2).Aicn; %chdata(i).sdd21f not scalled
        subplot(3,1,1)
        hold on
        plot(faxis_GHz, 20*log10(abs(PSXT*scale)),'r','Disp','PSXTK')
        icrxi=find(chdata(i).faxis >=param.fb/2,1,'first');
        subplot(3,1,2)
        grid on
        ILtemp=20*log10(abs(chdata(1).sdd21f));
        IL4ICR=interp1(chdata(1).faxis,ILtemp,chdata(i).faxis);
        scale=1/chdata(2).Aicn; % chdata(i).sdd21f not scalled
        ICR=-20*log10(abs(PSXT*scale))+IL4ICR;
        semilogx(faxis_GHz, ICR,'Disp', 'ICR')
        hold on
        stem(faxis_GHz(icrxi), ICR(icrxi),'g', 'disp', 'f_{Baud}/2')
    end
    subplot(3,1,1)
    title([param.base ' Losses']); ylabel('dB'); xlabel('GHz')
    grid on;  legend show
    subplot(3,1,2)
    title([param.base ' ICR']); ylabel('dB'); xlabel('GHz')
    ylim([0 80])
    xlim([.1 100])
    grid on; %legend show
    subplot(3,1,3)
    title([param.base ' ILD']); ylabel('dB'); xlabel('GHz')
    ylim([-3 3])
    grid on; legend show
end
% find loss values by interpolation using full data, no indexing - Adee 2022-08-28
total_loss = interp1(chdata(i).faxis, -20*log10(abs(chdata(1).sdd21)), 1/param.ui/2);
d2d_loss = interp1(chdata(i).faxis, -20*log10(abs(chdata(1).sdd21p_nodie)), 1/param.ui/2); 
VIP_VMP_loss = interp1(chdata(i).faxis, -20*log10(abs(chdata(1).sdd21p)), 1/param.ui/2); 
output_args.VTF_loss_dB_at_Fnq=total_loss;
output_args.IL_db_die_to_die_at_Fnq=d2d_loss;
output_args.VIP_to_VMP_IL_dB_at_Fnq=VIP_VMP_loss;
function [ V0 ] = FFE( C , cmx,spui, V )
% C      FFE taps
% cmx    number of precursors taps
% spui   samples per ui
% V      input signal
%speed ups implemented:
%1) ignore C(i)=0.  No need to circshift and then multiply by 0
%2) change ishift to shift an extra -cmx.  This avoids extra circshift at the end

V0=0;
if iscolumn(V); V=V.';end
for i=1:length(C)
    if C(i)~=0
        ishift=(i-1-cmx)*spui;
        V0=circshift(V',[ishift,0])*C(i)+V0;
    end
end
%V0=circshift(V0,[(-cmx)*spui,0]);
% disp(max(V0));


% begin yasuo patch 12/11/2018
% calculate sigma (standard deviation) value of PDF
function [ V0 ] = FFE_Fast( C,V_shift )
% C      FFE taps
% V      input signal separated into length(C) columns with circshift already performed 
% This function is only to speed up FFE in optimize_fom.  Since the signal that is being
% shifted is the same for all loops of TXFFE taps, a lot of time can be
% saved by pre-shifting it and remembering it across loops
% Another speed up:  only multiply by indices of C that are not 0

V0=0;
for i=1:length(C)
    if C(i)~=0
        V0=V_shift(:,i)*C(i)+V0;
    end
end
function [ V0 ] = Fract_T_FFE( V , skew_step)
% skew_step   sub UI skew assuming param.samples_per_ui
% V      input signal
% V0     output signal
% Richard Mellitz 8/17/2021
V0=0;
if iscolumn(V); V=V.';end
ishift=skew_step;
V0=circshift(V',[ishift,0])'+V;
V0=V0/2;
function out=Full_Grid_Matrix(in)

%create a full grid matrix of input variables
%used to create the full grid of all txffe cases
%example:
%Full_Grid_Matrix({ [1 2] [100 200] })
%out =
%     1   100
%     1   200
%     2   100
%     2   200
%
%input can also be mixed between numeric and cell of char
%example:
%Full_Grid_Matrix({ [1 2] {'A' 'B'} })
%out =
%    {[1]}    {'A'}
%    {[1]}    {'B'}
%    {[2]}    {'A'}
%    {[2]}    {'B'}

if ~iscell(in)
    error('input must be cell array of individual sweep variables');
end

num_columns=length(in);
num_cases=prod(cellfun('length',in));

cell_output=0;
cell_indices=cellfun(@(x) iscell(x),in);
if any(cell_indices)
    cell_output=1;
end
if cell_output
    for k=find(~cell_indices)
        in{k}=num2cell(in{k});
    end
end

if cell_output
    out=cell(num_cases,num_columns);
else
    out=zeros(num_cases,num_columns);
end

%num_repetitions controls how many times each element of the column
%repeats.  The first column is always just a copy of itself since every
%case will vary.
num_repetitions=1;
for k=num_columns:-1:1
    this_column=in{k}(:);
    %copy the column into a matrix to create the repetitions needed
    B=repmat(this_column,[1 num_repetitions]);
    %reshape into single column (actual repetitions)
    C=reshape(B',[numel(B) 1]);
    %repeat the single column to build the entire length required
    num_repeats=num_cases/length(C);
    D=repmat(C,[num_repeats 1]);
    out(:,k)=D;
    %determine how many repetitions the next column needs
    num_repetitions=num_repetitions*length(this_column);
end
function pdf=Init_PDF_Fast( EmptyPDF, values, probs)
%  p=cpdf(type, ...)
%
% CPDF is a probability mass function for discrete distributions or an
% approxmation of a PDF for continuous distributions.
%
% cpdf is internally normalized so that the sum of probabilities is 1
% (regardless of bin size).

% Internal fields:
% Min: *bin number* of minimum value.
% BinSize: size of PDF bins. Bin center is the representative value.
% Vec: vector of probabilities per bin.

pdf=EmptyPDF;

rounded_values_div_binsize=round(values/pdf.BinSize);
%values=pdf.BinSize*rounded_values_div_binsize;

% %speed up for small values round to 0 (because they are all much smaller than binsize)
% if all(values==0)
%     return;
% end
% 
% %speed up for all values rounded to the same bin
% %The output pdf is the same as the
% %empty pdf, but the x value is non-zero (but still scalar)
% if all(values==values(1))
%     pdf.Min=rounded_values_div_binsize(1);
%     pdf.x=values(1);
%     return;
% end
% 
% %The code below requires that values is
% %sorted.  Generally this should be true, but check to be sure
% if ~issorted(values)
%     [values,si]=sort(values);
%     rounded_values_div_binsize=rounded_values_div_binsize(si);
%     probs=probs(si);
% end


%pdf.x=values(1):pdf.BinSize:values(end);
pdf.x=pdf.BinSize*rounded_values_div_binsize(1):pdf.BinSize:pdf.BinSize*rounded_values_div_binsize(end);
pdf.Min=rounded_values_div_binsize(1);

pdf.y=zeros(size(pdf.x));
%The rounded values divided by binsize will reveal the bin number if
%pdf.Min is subtracted from it
bin_placement=rounded_values_div_binsize-pdf.Min+1;
%Can avoid one addition by inserting the first probability
%actually helps when calling this 2 million times
pdf.y(bin_placement(1))=probs(1);
for k=2:length(values)
    pdf.y(bin_placement(k)) = pdf.y(bin_placement(k))+probs(k);
end


%Have already ensured that sum(pdf.y)=1
%pdf.y=pdf.y/sum(pdf.y);

% if any(~isreal(pdf.y)) || any(pdf.y<0)
%     error('PDF must be real and nonnegative');
% end

% pMax=pdf.Min+length(pdf.y)-1;
% pdf.x = values(1):pdf.BinSize:pMax*pdf.BinSize;
function [MLSE_results] = MLSE(param,alpha,A_s,A_ni,PDF,CDF)
% OP.MLSE= 1 ... COM and VEC will be adjusted with MLSE CDF  
% OP.MLSE= 2 ... COM and VEC will be adjusted with MLSE Gaussian assumptions  
% Based on oif2022.580.00 / IEEE802. shakiba_3dj_01_230116 by Hossein Shakiba

qfuncinv = @(x) sqrt(2)*erfcinv(2*x);
qfunc = @(x) 0.5*erfc(x/sqrt(2));

%% step 0
    COM_from_matlab=20*log10(A_s/A_ni);
    L=param.levels;
    DER0=param.specBER;
%% step 1 from slide 6/5
    A_peak=(L-1)*A_s; % slide 6 A_s is main in appendix a
    main=A_peak;
    k_DER=qfuncinv(param.specBER);
    sigma_noise=sqrt(sum(PDF.y.*PDF.x.^2));
    SNR_dB=10*log10( 1/3*(L+1)/(L-1)*(A_peak^2)/sigma_noise^2) ;
    COM=SNR_dB-10*log10((L^2-1)/3*k_DER^2);
%     sprintf('COM from Matlab %g dB\n COM from slide 6 using Gaussian asumptions %g dB\n', COM_from_matlab ,COM)
if A_s >= A_ni
%% step 2 slide 10/8
    SNR_DFE=1/3*(L+1)/(L-1)*(A_peak^2)/sigma_noise^2;
%% step 2 slide 10/8
%     DER_DFE=    2/ ( L/(L-1) -qfunc(   (1-2*alpha)*main/(L-1)/sigma_noise )   )*(qfunc(main/(L-1)/sigma_noise));
%     DER_DFE_CDF=2/ ( L/(L-1)-CDF_ev(   (1-2*alpha)*main/(L-1),PDF,CDF )       )*CDF_ev((main/(L-1)),PDF,CDF);
%% step 3 side 11/9
    j=1:200;
    DER_MLSE=2*sum( j .* ((L-1)/L).^j .* qfunc( sqrt(1+(j-1)*(1-alpha)^2+alpha^2).* main/((L-1)*sigma_noise )       ));
    DER_MLSE_CDF=0; jj=1;
    DER_delta = inf;
    while DER_delta > .001
        last_DER_MLSE_CDF=DER_MLSE_CDF;
        DER_MLSE_CDF=2*( jj .* ((L-1)/L).^jj .* CDF_ev( sqrt(1+(jj-1)*(1-alpha)^2+alpha^2).* main/((L-1) ),PDF,CDF       ))+DER_MLSE_CDF;
        DER_delta= 1-last_DER_MLSE_CDF/DER_MLSE_CDF;
        jj=jj+1;
    end
%% step 4 slide 12/10
    SNR_DFE_eqivalent=SNR_DFE*(...
        (L-1)*sigma_noise/main *    qfuncinv(...
        1/2 *DER_MLSE*(L/(L-1) - qfunc((1-2*alpha)*main/(L-1)*sigma_noise )) ...
        )   ...
        )^2;
    SNR_DFE_eqivalent_CDF=SNR_DFE*(...
        (L-1)/main *    CDF_inv_ev(...
        1/2 *DER_MLSE_CDF*(L/(L-1) - CDF_ev((1-2*alpha)*main/(L-1),PDF,CDF )) ...
        ,PDF, CDF )   ...
        )^2;

%% step 5 slide 13/11
delta_com=10*log10(SNR_DFE_eqivalent/SNR_DFE);
delta_com_CDF=10*log10(SNR_DFE_eqivalent_CDF/SNR_DFE);
new_com_CDF=COM_from_matlab+delta_com_CDF;
else
    warning('MLSE not applied because there is more noise than signal')
    DER_MLSE=[];
    DER_MLSE_CDF=[];
    SNR_DFE_eqivalent=[];
    SNR_DFE_eqivalent_CDF=[];
    new_com_CDF=COM_from_matlab;
    delta_com_CDF=0;
    delta_com=0;
    SNR_DFE=[];
end

%%
MLSE_results.COM_from_matlab=COM_from_matlab;
MLSE_results.SNR_DFE=SNR_DFE;
MLSE_results.DER_MLSE_Gaussian=DER_MLSE;   
MLSE_results.DER_MLSE_CDF=DER_MLSE_CDF;
MLSE_results.sigma_noise=sigma_noise;
MLSE_results.SNR_dB=SNR_dB ;
MLSE_results.SNR_DFE_eqivalent_Gaussian=SNR_DFE_eqivalent;
MLSE_results.SNR_DFE_eqivalent_CDF=SNR_DFE_eqivalent_CDF;
MLSE_results.COM_Gaussian=new_com_CDF;
MLSE_results.COM_CDF=new_com_CDF;
MLSE_results.k_DER=k_DER;
MLSE_results.delta_com_CDF=delta_com_CDF;
MLSE_results.delta_com_Gaussian=delta_com;



function [MLSE_results] = MLSE_instu(param,alpha,A_s,A_ni,PDF,CDF)
% OP.MLSE= 1 ... COM and VEC will be adjusted with MLSE CDF
% OP.MLSE= 2 ... COM and VEC will be adjusted with MLSE Gaussian assumptions
% Based on oif2022.580.00 / IEEE802. shakiba_3dj_01_230116 by Hossein Shakiba

qfuncinv = @(x) sqrt(2)*erfcinv(2*x);
qfunc = @(x) 0.5*erfc(x/sqrt(2));

%% step 0
COM_from_matlab=20*log10(A_s/A_ni);
L=param.levels;
DER0=param.specBER;
%% step 1 from slide 6/5
A_peak=(L-1)*A_s; % slide 6 A_s is main in appendix a
main=A_peak;
k_DER=qfuncinv(param.specBER);
sigma_noise=sqrt(sum(PDF.y.*PDF.x.^2));
SNR_dB=10*log10( 1/3*(L+1)/(L-1)*(A_peak^2)/sigma_noise^2) ;
COM=SNR_dB-10*log10((L^2-1)/3*k_DER^2);
%     sprintf('COM from Matlab %g dB\n COM from slide 6 using Gaussian asumptions %g dB\n', COM_from_matlab ,COM)
if A_s >= A_ni
    %% step 2 slide 10/8
    SNR_DFE=1/3*(L+1)/(L-1)*(A_peak^2)/sigma_noise^2;
    %% step 2 slide 10/8
    snr_dfe = @(der,PDF,CDF) -10*log10((A_s./CDF_inv_ev(der,PDF,CDF)).^2)+10*log10((L^2-1)/3*qfuncinv(der).^2) ;

    %% step 3 side 11/9
    j=1:200;
    DER_MLSE=2*sum( j .* ((L-1)/L).^j .* qfunc( sqrt(1+(j-1)*(1-alpha)^2+alpha^2).* main/((L-1)*sigma_noise )       ));
    DER_MLSE_CDF=0; jj=1;
    DER_delta = inf;
    while DER_delta > .001
        last_DER_MLSE_CDF=DER_MLSE_CDF;
        DER_MLSE_CDF=2*( jj .* ((L-1)/L).^jj .* CDF_ev( sqrt(1+(jj-1)*(1-alpha)^2+alpha^2).* main/((L-1) ),PDF,CDF       ))+DER_MLSE_CDF;  
        DER_delta= 1-last_DER_MLSE_CDF/DER_MLSE_CDF;
        jj=jj+1;
    end
    %%
    dscale=.05;
    scale=1;
    last_scale_tune=inf;
    scale_tune=inf;
    while abs(scale_tune) >= .1
        istart=-PDF.Min+1;
        scale=scale-dscale;
        PDF_SCALED = scalePDF(PDF,scale);
        cdf_scaled=pdf_to_cdf(PDF_SCALED);
        test_snr=snr_dfe(DER_MLSE_CDF,PDF_SCALED,cdf_scaled.y);
        test_DER=2*(L-1)/L* (CDF_ev(main/(L-1),PDF_SCALED,cdf_scaled.y) );
        scale_tune=(test_DER-DER_MLSE_CDF)/DER_MLSE_CDF;
        if sign(scale_tune) ~= sign(last_scale_tune)
            % scale=scale+dscale % back up
            dscale=-dscale/2;
        end
        last_scale_tune=scale_tune;
    end
    new_com_CDF=10*log10((A_s./CDF_inv_ev(DER0,PDF_SCALED,cdf_scaled.y)).^2);
    delta_com=new_com_CDF-10*log10((A_s./CDF_inv_ev(DER0,PDF,CDF)).^2); 
else
    warning('MLSE not applied because there is more noise than signal')
    DER_MLSE=[];
    DER_MLSE_CDF=[];
    SNR_DFE_eqivalent=[];
    SNR_DFE_eqivalent_CDF=[];
    new_com_CDF=COM_from_matlab;
    delta_com_CDF=0;
    delta_com=0;
    SNR_DFE=[];
    PDF_SCALED=[];
    cdf_scaled=[];
end

%%
MLSE_results.COM_from_matlab=COM_from_matlab;
MLSE_results.DER_MLSE_Gaussian=DER_MLSE;
MLSE_results.DER_MLSE_CDF=DER_MLSE_CDF;
MLSE_results.sigma_noise=sigma_noise;
MLSE_results.k_DER=k_DER;
MLSE_results.COM_CDF=new_com_CDF;
MLSE_results.delta_com_CDF=delta_com;
MLSE_results.delta_com_Gaussian=delta_com;
MLSE_results.PDF=PDF_SCALED;
MLSE_results.CDF=cdf_scaled.y;
MLSE_results.PDF_scale=scale;




function MMSE_results = MMSE(PSD_results,sbr, cursor_i ,param, OP ) ;
if 1
    num_ui=param.num_ui_RXFF_noise;
    M=param.samples_per_ui;
    L=param.levels;
    sigma_X2=(L^2-1)/(3*(L-1)^2);
    fb=param.fb;
    R_LM=param.R_LM;
end
h=sbr(mod(cursor_i,M)+1:end-mod(cursor_i,M)); % align to sample point
h=reshape(h,1,[]); % make row vectors
h=[ h(1:floor(length(h)/M)*M) ];
h= [h zeros(1,num_ui*M-length(h)) ];
h=h(1:M:end);% resample
N=length(h);
dw=param.RxFFE_cmx ; % equalizer precuror tapsindx(1:N)=(1:N)-5-1;
dh=(cursor_i-mod(cursor_i,M))/M ; % precuror taps in h
if param.N_bg == 0
    Nw= param.RxFFE_cmx+1+param.RxFFE_cpx; %  total number of equalizer taps
    bmax=param.bmax;
    bmin=param.bmin ;
    wmax= [   ones(1,param.RxFFE_cmx-1)*param.ffe_tapn_max   param.ffe_pre_tap1_max 1.0  param.ffe_post_tap1_max     ones(1,param.RxFFE_cpx-1)*param.ffe_tapn_max  ];
    wmin= [  -ones(1,param.RxFFE_cmx-1)*param.ffe_tapn_max  -param.ffe_pre_tap1_max 1.0  -param.ffe_post_tap1_max   -ones(1,param.RxFFE_cpx-1)*param.ffe_tapn_max  ];
    idx=[];
else
    Nfloating_taps=param.N_bf*param.N_bg;
    Nmax=param.N_bmax;
    Nfix=param.RxFFE_cmx+1+param.RxFFE_cpx;
    Ng=param.N_bg;
    Nf=param.N_bf;
    Nw= dw+Nmax+1;%  total span of equalizer taps including floating taps
    Nwft=param.RxFFE_cmx+1+param.RxFFE_cpx+Nfloating_taps;%  total number of equalizer taps including floating taps
    hisi=h(dh+2:((dh-dw)+Nw));
    [idx]=findbankloc( hisi ,param.N_tail_start,param.N_bmax,Nf,inf,param.bmaxg,Ng  ); % using maximum power in hisi
    idx=sort(idx);
    bmax=param.bmax;
    bmin=param.bmin ;
    wmax= [   ones(1,param.RxFFE_cmx-1)*param.ffe_tapn_max   param.ffe_pre_tap1_max 1.0  param.ffe_post_tap1_max   ones(1,param.RxFFE_cpx-1)*param.ffe_tapn_max  ones(1,Nfloating_taps)*param.bmaxg ];
    wmin= [  -ones(1,param.RxFFE_cmx-1)*param.ffe_tapn_max  -param.ffe_pre_tap1_max 1.0  -param.ffe_post_tap1_max -ones(1,param.RxFFE_cpx-1)*param.ffe_tapn_max -ones(1,Nfloating_taps)*param.bmaxg  ];
    if 0
        figure
        set(gcf, 'tag', 'COM');%movegui(gcf,'south');
        hindx=1:length(hisi);
        stem(hindx,hisi)
        hold on
        stem(idx,hisi(idx))
    end

end
Nb=param.ndfe; % DFE taps
d=dw+dh; % used for index in algorithms
indx(1:N)=(1:N)-dh-1;

S_n=PSD_results.S_n; % total agregate noise PSD

Rn=ifft(S_n)*fb;
%% HH and R
Rnn=toeplitz(Rn(1:Nw),Rn(1:Nw));
hc1=[ h  zeros(1,Nw-1) ];
hr1=[ h(1) zeros(1,Nw-1)];
H=toeplitz(hc1,hr1);
if param.N_bg ~= 0
    H=H( :,[1:Nfix idx]);
end
HH= H'*H;
if param.N_bg ~= 0
 Rnn=Rnn( [1:Nfix idx-param.N_tail_start+1+Nfix],[1:Nfix idx-param.N_tail_start+1+Nfix]);
end
R=HH+Rnn/sigma_X2;
%% hb and h0
Hb= H(d+2:d+Nb+1,:);
h0=H(d+1,:);
% display(floor(h0));

%% Ib and zb (slide 10)
ib=eye(Nb);
zb=zeros(1,Nb);
wbl= [ R  -Hb' -h0';...
    -Hb  ib  zb'; ...
    h0  zb   0]\[h0'; zb' ;1];

%% re-adjust Nw to number of used taps
if param.N_bg ~= 0
    Nw=Nwft;
end
%% check equalized pulse
w=wbl(1:Nw);
b=wbl(Nw+1:length(wbl)-1); % dfe taps before limits are applied

%% apply blim (slide 11)  <---- need help here How do I get to RxFFE tap coefficents, C?


blim = min(bmax(:), max(bmin(:), b));
if (Nb > 0) && ~isequal(b, blim)
    wl = [R, -h0'; h0, 0]\[h0'+Hb'*blim; 1];
    w = wl(1:Nw);
end

wlim = min(wmax(:)*w(1+dw), max(wmin(:)*w(1+dw), w));
if ~isequal(w, wlim)
    wlim = wlim/(h0*wlim); % Ensure the equalized pulse amplitude is 1.
    if Nb > 0
        b = Hb*wlim; % Update the feedback coefficients.
        blim = min(bmax(:), max(bmin(:), b));
    end
    wl = [R, -h0'; h0, 0]\[h0'+Hb'*blim; 1];
    wl = wl(1:Nw);
    w = min(wmax(:)*wl(1+dw), max(wmin(:)*wl(1+dw), wl));
end

w=w(1:Nw) ;
% w=wl(1:Nw); % RxFFE taps before limits are applied
Craw=w/w(dw+1); % returned Rx FFE taps
% re-align Cmod to floating tap locations
if param.N_bg ~= 0
    C=Craw;
    C(Nfix+1:Nmax)=0;
    C(idx-param.N_tail_start+1+Nfix)=Craw(Nfix+(1:Nfloating_taps));
else
    C=Craw;
end
if 0
    figure
    set(gcf, 'tag', 'COM');movegui(gcf,'south');
     subplot(3,1,1)
     stem(-dw:Nmax-dw-1,C,'DisplayName','Returned RxFFE taps')
     title('solved Rxffe C')
    legend show
    subplot(3,1,2)
    resp=conv(C,h);
    resp=resp(dw+1:N+dw);
    stem(indx(1:length(resp)),resp,'disp', 'Rxffe filtered h with C')
    title('convolved response with C')
    legend show
    subplot(3,1,3)
    resp=conv(w,h);
    resp=resp(dw+1:N+dw);
    stem(indx(1:length(resp)),resp,'disp', 'Rxffe filtered h with w(ts) not normalized')
    title('convolved response with  w')
    legend show
end

MMSE_results.floating_tap_locations=idx;
MMSE_results.C=C;
MMSE_results.sigma_e=(w'*R*w+1+b'*b-2*w'*h0'-2*w'*Hb'*b);
MMSE_results.FOM=20*log10((R_LM/(L-1)/MMSE_results.sigma_e));



function output_args=Output_Arg_Fill(output_args,sigma_bn,Noise_Struct,COM_SNR_Struct,param,chdata,fom_result,OP)

%not all output_args are filled here but most are

switch lower(OP.TDECQ)
    case { false 'none' } % should be the default
        output_args.VMA=[];
    case 'vma'
        est_vma=vma(fom_result.sbr,param.samples_per_ui);
        output_args.VMA=est_vma.VMA;
    otherwise
        error('%s not recognized for Histogram_Window_Weigh this feature is limited',OP.TDECQ)
end

fileset_str=str2csv({chdata.base});
output_args.file_names=sprintf('"%s"', fileset_str);
% [ahealey] Echo the termination parameters in the output arguments..
for odt_param = {'R_diepad', 'C_diepad', 'L_comp', 'C_bump'}
    output_args.(odt_param{:}) = param.(odt_param{:});
end
% [ahealey] End of modifications.
for pkg_params = {'levels', 'Pkg_len_TX', 'Pkg_len_NEXT', 'Pkg_len_FEXT', 'Pkg_len_RX','R_diepad','pkg_Z_c','C_v'}
    output_args.(pkg_params{:})= param.(pkg_params{:});
end
output_args.baud_rate_GHz=param.fb/1e9;
output_args.f_Nyquist_GHz = param.fb/2e9;
output_args.BER=param.specBER;
output_args.FOM = fom_result.FOM;
output_args.sigma_N=Noise_Struct.sigma_N;
output_args.DFE4_RSS=norm(fom_result.DFE_taps(4:end));
output_args.DFE2_RSS=norm(fom_result.DFE_taps(2:end));
output_args.tail_RSS=fom_result.tail_RSS;
output_args.channel_operating_margin_dB=COM_SNR_Struct.COM;
output_args.available_signal_after_eq_mV=1000*COM_SNR_Struct.A_s;
output_args.peak_uneq_pulse_mV=1000*max(abs(chdata(1).uneq_pulse_response));
try
    output_args.uneq_FIR_peak_time=chdata(1).t(chdata(1).uneq_imp_response==max(chdata(1).uneq_imp_response));
catch
    output_args.uneq_FIR_peak_time=[];
end
output_args.steady_state_voltage_mV = 1000*fom_result.A_f; % RIM 7/03/2019 use peak from optimize_FOM
its=find(chdata(1).eq_pulse_response>=max(chdata(1).eq_pulse_response),1,'first');
isumend=min(its+param.N_v*param.samples_per_ui,length(chdata(1).eq_pulse_response));
output_args.steady_state_voltage_weq_mV = 1000*sum(chdata(1).eq_pulse_response(1:isumend) )/param.samples_per_ui;

if OP.RX_CALIBRATION== 1
    output_args.sigma_bn=sigma_bn;
else
    output_args.sigma_bn=[];
end
output_args.Peak_ISI_XTK_and_Noise_interference_at_BER_mV=1000*COM_SNR_Struct.A_ni;
output_args.peak_ISI_XTK_interference_at_BER_mV=1000*Noise_Struct.peak_interference_at_BER;
output_args.peak_ISI_interference_at_BER_mV=1000*Noise_Struct.thru_peak_interference_at_BER;
output_args.equivalent_ICI_sigma_assuming_PDF_is_Gaussian_mV=Noise_Struct.sci_sigma*1000;

if OP.RX_CALIBRATION == 0
    output_args.peak_MDXTK_interference_at_BER_mV=1000*Noise_Struct.crosstalk_peak_interference_at_BER;
    output_args.peak_MDNEXT_interference_at_BER_mV=1000*Noise_Struct.MDNEXT_peak_interference;
    output_args.peak_MDFEXT_interference_at_BER_mV=1000*Noise_Struct.MDFEXT_peak_interference;
else
    output_args.peak_MDXTK_interference_at_BER_mV=[];
    output_args.peak_MDNEXT_interference_at_BER_mV=[];
    output_args.peak_MDFEXT_interference_at_BER_mV=[];
end
%output_args.ICN_mV=ICN*1000;
%             output_args.ICN_test_mV=ICN_test*1000;
xtk=param.num_next+param.num_fext;
if xtk>0 && OP.RX_CALIBRATION ==0 && OP.TDMODE==0
    %output_args.MDNEXT_ICN_92_46_mV=MDNEXT_ICN*1000;
    %output_args.MDFEXT_ICN_92_47_mV=MDFEXT_ICN*1000;
    output_args.equivalent_ICN_assuming_Gaussian_PDF_mV=Noise_Struct.cci_sigma*1000;
else
    output_args.MDNEXT_ICN_92_46_mV=0;
    output_args.MDFEXT_ICN_92_47_mV=0;
    output_args.equivalent_ICN_assuming_PDF_is_Gaussian_mV=0;
end
%output_args.SNR_ISI_XTK_normalized_1_sigma=20*log10(A_s/(peak_interference_at_BER/qfuncinv(param.specBER))); modified by Yasuo Hidaka, 8/7/17
if 1
    output_args.SNR_ISI_XTK_normalized_1_sigma=20*log10(COM_SNR_Struct.A_s/(Noise_Struct.peak_interference_at_BER/sqrt(2)/erfcinv(2*param.specBER)));
    output_args.SNR_ISI_est=fom_result.SNR_ISI;
    output_args.Pmax_by_Vf_est=fom_result.Pmax_by_Vf;
    output_args.Tr_measured_from_step_ps=fom_result.Tr_measured_from_step/1e-12;
end


switch param.CTLE_type
    case 'CL93'
        output_args.CTLE_zero_poles=[param.CTLE_fz(fom_result.ctle) param.CTLE_fp2(fom_result.ctle) param.CTLE_fp1(fom_result.ctle)];
        output_args.CTLE_DC_gain_dB=param.ctle_gdc_values(fom_result.ctle);
        output_args.g_DC_HP=[];
        output_args.HP_poles_zero=[];
    case 'CL120d'
        output_args.CTLE_zero_poles=[param.CTLE_fz(fom_result.ctle) param.CTLE_fp2(fom_result.ctle) param.CTLE_fp1(fom_result.ctle)];
        output_args.CTLE_DC_gain_dB=param.ctle_gdc_values(fom_result.ctle);
        output_args.g_DC_HP=param.g_DC_HP_values(fom_result.best_G_high_pass);
        output_args.HP_poles_zero=param.f_HP(fom_result.best_G_high_pass);
    case 'CL120e'
        output_args.CTLE_zero_poles=[param.CTLE_fz(fom_result.ctle) param.f_HP_Z(fom_result.ctle) param.CTLE_fp2(fom_result.ctle) param.CTLE_fp1(fom_result.ctle) param.f_HP_P(fom_result.ctle)];
        output_args.CTLE_DC_gain_dB=param.ctle_gdc_values(fom_result.ctle);
        output_args.g_DC_HP=[];
        output_args.HP_poles_zero=[];
end
output_args.TXLE_taps=fom_result.txffe;
if  length(output_args.TXLE_taps) >= 3
    output_args.Pre2Pmax = -output_args.TXLE_taps(end-2)/output_args.TXLE_taps(end-1);
else
    output_args.Pre2Pmax=[];
end
output_args.DFE_taps=fom_result.DFE_taps;
if param.Floating_DFE ||  param.Floating_RXFFE
    output_args.floating_tap_locations=fom_result.floating_tap_locations;
else
    output_args.floating_tap_locations=[];
end

if OP.RxFFE
    output_args.RxFFE=fom_result.RxFFE;
    output_args.RxFFEgain=param.current_ffegain;
else % Yasou Hidaka 11/20/2018 help to align csv file columns
    output_args.RxFFE=[];
    output_args.RxFFEgain=[];
end

output_args.itick=fom_result.itick;

% Calculation of error propagation and burst probability
if OP.nburst>0
    [p_burst,p_error_propagation]=Burst_Probability_Calc(COM_SNR_Struct,fom_result.DFE_taps,param,OP);
    output_args.error_propagation_probability = p_error_propagation;
    output_args.burst_probabilities = p_burst;
else
    output_args.error_propagation_probability = [];
    output_args.burst_probabilities = [];
end


%begin yasuo patch 12/11/2018
% collect sigma values to report
% pdf2sgm() is a function to calculate sigma value from PDF
% It is added at the end of this file code.
% I am not sure if an equivalent function already exists.
output_args.sgm_Ani__isi_xt_noise = pdf2sgm(COM_SNR_Struct.combined_interference_and_noise_pdf);
output_args.sgm_isi_xt = pdf2sgm(Noise_Struct.isi_and_xtalk_pdf);
output_args.sgm_noise__gaussian_noise_p_DD  = pdf2sgm(Noise_Struct.noise_pdf);
output_args.sgm_p_DD = pdf2sgm(Noise_Struct.p_DD);
output_args.sgm_gaussian_noise = pdf2sgm(Noise_Struct.gaussian_noise_pdf);
output_args.sgm_G = Noise_Struct.sigma_G;
output_args.sgm_rjit = Noise_Struct.sigma_rjit;
output_args.sgm_N = Noise_Struct.sigma_N;
output_args.sgm_TX = Noise_Struct.sigma_TX;
output_args.sgm_isi = pdf2sgm(Noise_Struct.sci_pdf);
if OP.RX_CALIBRATION == 0
    output_args.sgm_xt  = pdf2sgm(Noise_Struct.cci_pdf);
else
    output_args.sgm_xt=[];
end
% end yasuo patch

% r259 putting COM, VEO and loss last in report
%         output_args.VEO_normalized = (A_s-A_ni)/A_s;
output_args.VEC_dB = COM_SNR_Struct.VEC_dB;
output_args.VEO_mV = COM_SNR_Struct.VEO_mV;
if OP.RX_CALIBRATION ==0 && OP.EW == 1
    output_args.EW_UI_est=COM_SNR_Struct.EW_UI;
    output_args.eye_contour=COM_SNR_Struct.eye_contour;
    output_args.VEO_window_mUI= param.T_O;
else
    output_args.EW_UI_est=[];
    output_args.eye_contour=[];
    output_args.VEO_window_mUI= [];
end

if sum(param.AC_CM_RMS) ~= 0
    output_args.sigma_ACCM_at_tp0_mV=chdata(1).sigma_ACCM_at_tp0*1000;
    fprintf(' AC RMS at TP0 = %.3g mV \n',output_args.sigma_ACCM_at_tp0_mV)
    output_args.sigma_AC_CCM_at_rxpkg_output_mV=chdata(1).CD_CM_RMS*1000; %
else
    output_args.sigma_ACCM_at_tp0_mV=[];
    output_args.sigma_AC_CCM_at_rxpkg_output_mV=[];
end
if OP.MLSE
    output_args.COM_orig=COM_SNR_Struct.COM_orig;
    output_args.VEC_dB_orig=COM_SNR_Struct.VEC_dB_orig;
end
%
output_args.COM_dB=COM_SNR_Struct.COM;
% end yasuo patch
% begin yasuo patch 3/18/2019
output_args.DER_thresh = COM_SNR_Struct.threshold_DER;
% end yasuo patch
function [ seq syms syms_nrz ] = PRBS13Q( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


taps = ([13 12 2  1]);
seed =([0 0 0 0 0 1 0 1 0 1 0 1 1]);
[seq_nrz c] =LFSR(seed,taps);
seq_nrz=2*(seq_nrz-0.5);
seq=pam(seq_nrz);
% syms=round(2*(seq+1));
syms((round(2*(seq+1))/2==2))=3;
syms((round(2*(seq+1))/2==1.5))=2;
syms((round(2*(seq+1))/2==.5))=1;
syms((round(2*(seq+1))/2==0))=0;

% syms_nrz=((seq_nrz+1)/2);

syms_nrz=seq_nrz;


function[seq c]=LFSR(s,t)
%s=initial state of LFSR, you can choose any lenght of LFSR
%Instruction:========== 
%Save LFSR.m in your current directory and type following
%on Command window for simulating 5 bit LFSR with tap [5 2]
%---------------------
%>>s=[1 1 0 0 1] 
%>>t=[5 2]
%>>[seq c] =LFSR(s,t) 
%---------------------------
%seq = generated sequence
%c will be matrix containing the states of LFSR raw wise
% 
%-----------------------------------------------------------
%If any doubt, confusion or feedback please contact me
% NIKESH BAJAJ 
% bajaj.nikkey@gmail.com (+91-9915522564)
% Asst. Professor at Lovely Profesional University
% Masters from Aligarh Muslim University,INDIA 
%--------------------------------------------------
n=length(s);
c(1,:)=s;
m=length(t);
for k=1:2^n-2;
b(1)=xor(s(t(1)), s(t(2)));
if m>2;
    for i=1:m-2;
    b(i+1)=xor(s(t(i+2)), b(i));
    end
end
j=1:n-1;
s(n+1-j)=s(n-j);
s(1)=b(m-1);
c(k+1,:)=s;
end
seq=c(:,n)';

function [ dataout ] = pam( data )
% mapping data usng Grey Coding
for i=1:2:floor(length(data)/2)*2
    if data(i:i+1)==[ -1 -1 ]
        dataout(ceil(i/2)) = -1;
    elseif data(i:i+1)==[ -1 1 ]
        dataout(ceil(i/2)) = -1/3;
    elseif data(i:i+1)==[ 1 1 ]
        dataout(ceil(i/2)) = 1/3;
    elseif data(i:i+1)==[ 1 -1 ]
        dataout(ceil(i/2)) = 1;
    end
end
function RILN_TD_struct=RILN_TD(sdd21,RIL,faxis_f2,OP,param,A_T)
db = @(x) 20*log10(abs(x));
disp('computing TD_RILN...')
sdd21=squeeze(sdd21);
if  iscolumn(sdd21)
    sdd21=sdd21.';
end
RIL=squeeze(RIL);
if  iscolumn(RIL)
    RIL=RIL.';
end
print_for_codereview=1;
if exist('OP','var')
    X=sinc(faxis_f2*param.ui)*param.ui*1e9;
    
    H_bt=Bessel_Thomson_Filter(param,faxis_f2,1);
    H_bw=Butterworth_Filter(param,faxis_f2,1);
    H_t = exp(-(pi*faxis_f2/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
    H_tw=Tukey_Window(faxis_f2,param);
    H_tw=ones(1,length(faxis_f2) );
    [RILN_TD_struct.REF.FIR, ...
        RILN_TD_struct.REF.t, ...
        RILN_TD_struct.REF.causality_correction_dB, ...
        RILN_TD_struct.REF.truncation_dB] = s21_to_impulse_DC(sdd21.*H_bw.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    RILN_TD_struct.REF.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.REF.FIR);
    [RILN_TD_struct.FIT.FIR, ...
        RILN_TD_struct.FIT.t, ...
        RILN_TD_struct.FIT.causality_correction_dB, ...
        RILN_TD_struct.FIT.truncation_dB] = s21_to_impulse_DC(RIL.*H_bw.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    RILN_TD_struct.FIT.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.FIT.FIR);   
    ipeak=find(RILN_TD_struct.REF.PR==max(RILN_TD_struct.REF.PR),1,'first');
    NrangeUI=1000;
    range_end=min(ipeak+param.samples_per_ui*NrangeUI,min(length(RILN_TD_struct.FIT.FIR), length(RILN_TD_struct.REF.FIR) ) ) ; 
    range=ipeak:range_end;
    RILN_TD_struct.ILN=RILN_TD_struct.FIT.PR(range)-RILN_TD_struct.REF.PR(range);
    RILN_TD_struct.t=RILN_TD_struct.FIT.t(range);
    RILN_TD_struct.FOM=-inf;
    RILN_TD_struct.FOM_PDF=-inf;
    rms_fom=-inf;
    for im=1:param.samples_per_ui
        RILN_TD_struct.FOM=max(RILN_TD_struct.FOM, norm( RILN_TD_struct.ILN(im:param.samples_per_ui:end)));
        [ pdf ] = get_pdf_from_sampled_signal(  RILN_TD_struct.ILN(im:param.samples_per_ui:end), param.levels, OP.BinSize ,0);
        rms=sqrt(pdf.y*pdf.x(:).^2)*sqrt(2);
        cdf=pdf; cdf.y=cumsum(pdf.y);
        %         cursors = d_cpdf(OP.BinSize,param.a_thru*[-1 1], [1 1]/2);
        %         signal_and_isi_pdf = conv_fct(cursors, pdf);
        %         cdf=signal_and_isi_pdf; cdf.y=cumsum(signal_and_isi_pdf.y);
        if print_for_codereview % remove once all checked out
            h=figure(191);set(gcf,'Tag','COM');
            semilogy(-cdf.x,cdf.y);
%             xlim ([0,-cdf.x(1)])
            ylim([param.specBER 1]);title ('CDF of RILN')
            hold on
        end
        if rms>rms_fom
            rms_fom=rms;
            RILN_TD_struct.FOM_PDF= -cdf.x(find(cdf.y >= param.specBER, 1, 'first'));
            RILN_TD_struct.PDF=pdf;
        end
    end
    pdf_from_norm=normal_dist(RILN_TD_struct.FOM, 7 , OP.BinSize);
    RILN_TD_struct.SNR_ISI_FOM=db(RILN_TD_struct.FIT.PR(ipeak)/RILN_TD_struct.FOM);
    RILN_TD_struct.SNR_ISI_FOM_PDF=db(RILN_TD_struct.FIT.PR(ipeak)/RILN_TD_struct.FOM_PDF);
    if print_for_codereview % remove once all checked out
        figure(9003);set(gcf,'Tag','COM');
        plot(RILN_TD_struct.t,RILN_TD_struct.ILN,'disp','td riln')
        hold on
        plot(RILN_TD_struct.FIT.t,RILN_TD_struct.FIT.PR,'disp','fit')
        plot(RILN_TD_struct.REF.t,RILN_TD_struct.REF.PR,'disp','ref')
        hold off
        fprintf('SNR ISI FOM rms = %g dB;   SNR ISI FOM PDF = %g dB\n',RILN_TD_struct.SNR_ISI_FOM,RILN_TD_struct.SNR_ISI_FOM_PDF)
        figure(9004);set(gcf,'Tag','COM');
        semilogy(RILN_TD_struct.PDF.x,RILN_TD_struct.PDF.y,'disp','actual PDF')
        hold on
        semilogy(pdf_from_norm.x,pdf_from_norm.y,'disp','PDF using Gaussian assumed PDF');
        ylim([param.specBER max([RILN_TD_struct.PDF.y pdf_from_norm.y])]);title ('Compare actual PDF to Gaussian')
        grid on
        legend('show')
    end
end
function is_illegal=RXFFE_Illegal(C,param,last_index)

%check if RXFFE taps are illegal
%C = RXFFE taps
%param = COM param struct
%last_index is used when computing illegality prior to Backoff.  It will be set so taps
% in the Backoff region are not considered in the legality check.

%If last index is omitted, set it to length(C)
if nargin<3
    last_index=length(C);
end

is_illegal=0;

%Check cursor tap
Ccur_i=param.RxFFE_cmx+1;
if C(Ccur_i) < param.ffe_main_cursor_min
    is_illegal=1;
    return;
end

%Check postcursors
if param.ffe_post_tap_len ~=0
    if abs(C(Ccur_i +1)) > param.ffe_post_tap1_max
        is_illegal=1;
        return;
    end
    if (param.ffe_post_tap_len > 1)
        if sum(abs(C((Ccur_i +2):last_index))  > param.ffe_tapn_max)
            is_illegal=1;
            return;
        end
    end
end

%Check precursors
if param.ffe_pre_tap_len ~=0
    if abs(C(Ccur_i -1)) > param.ffe_pre_tap1_max
        is_illegal=1;
        return;
    end
    if (param.ffe_pre_tap_len > 1)
        %                                             if sum(abs(C((Ccur_i +2):end)) > param.ffe_tapn_max) , continue; end
        if sum(abs(C(1:(Ccur_i - 2))) > param.ffe_tapn_max)
            is_illegal=1;
            return;
        end % 11.22.2018 Yasou Hadaka
    end
end
function S =R_series2(zref,f,R)
r=ones(1,length(f))*R;
S.Parameters(1,1,:) =  r./(r + 2*zref);
S.Parameters(2,2,:) =  r./(r + 2*zref);
S.Parameters(2,1,:) = (2*zref)./(r + 2*zref);
S.Parameters(1,2,:) =  (2*zref)./(r + 2*zref);
% Sm=sparameters(S.Parameters,f,zref);

function H_tw=Raised_Cosine_Filter(param,f,use_RC)

if use_RC
    H_tw = Tukey_Window(f,param ,param.RC_Start, param.RC_end);% add tw filter;
else
    H_tw=ones(1,length(f));
end
function SLD=SL(S,f,R)
% source load impact return loss add to S21
% S and SLD are the same structure
% S.Parameters
% S.Impedance
% S.NumPorts
% S.Frequencies
SLD=S; % assign the fields
zref=100;
if R==0
    warndlg('Termination should not be set to zero');
    SLD=S;
    return
end

if R > zref
    spr =R_series2(zref,f,(R-zref)); % make series sparameter
    %     SLD=sparameters(cascadesparams(S.Parameters,spr.Parameters),f,S.Impedance); % Casdeade
    [ SLD.Parameters(1,1,:), SLD.Parameters(1,2,:), SLD.Parameters(2,1,:), SLD.Parameters(2,2,:)]  = ...
        combines4p( ...
        spr.Parameters(1,1,:), spr.Parameters(1,2,:), spr.Parameters(2,1,:), spr.Parameters(2,2,:),...
        S.Parameters(1,1,:),   S.Parameters(1,2,:),   S.Parameters(2,1,:),   S.Parameters(2,2,:) ...
        );
elseif R < zref
    spr =r_parrelell2(zref,f,-R*zref/(R-zref)); % make parrellel sparameter
    %     SLD=sparameters(cascadesparams(S.Parameters,spr.Parameters),f,S.Impedance);
    [ SLD.Parameters(1,1,:), SLD.Parameters(1,2,:), SLD.Parameters(2,1,:), SLD.Parameters(2,2,:)]  = ...
        combines4p( ...
        spr.Parameters(1,1,:),  spr.Parameters(1,2,:), spr.Parameters(2,1,:), spr.Parameters(2,2,:),...
        S.Parameters(1,1,:),   S.Parameters(1,2,:),   S.Parameters(2,1,:),   S.Parameters(2,2,:) ...
        );
else
    SLD=S;
end

%%

function S_RN_of_f = S_RN(f,G_DC,G_DC2,param)
p1=param.CTLE_fp1(1);
z1=param.CTLE_fz(1);
p2=param.CTLE_fp2(1);
zlf=param.f_HP(1);
plf=param.f_HP(1);
f_b=param.fb;
f_r=param.f_r;
eta_0=param.eta_0;
H_CTF =   ( 10^(G_DC/20)+ 1i*f/z1 ) .*( 10^(G_DC2/20) + 1i*f/zlf )./ ( ( 1+1i*f/p1) .* ( 1+1i*f/p2)  .* ( 1+1i*f/plf));
H_R =1./polyval([1 2.613126 3.414214 2.613126 1], 1i*f./(f_r*f_b));
S_RN_of_f = eta_0/2.*abs( H_CTF.*H_R).^2; %EQ healey_3dj_01_2401 slide 5
if 0
    figure
    set(gcf, 'tag', 'COM');movegui(gcf,'southeast');
    % see if it looks correct
    semilogx(f/1e9, 20*log10( abs( H_CTF.*H_R)) );
    ylabel('dB');
    xlabel('GHz');
    title( 'H_ctf with H_r')
    grid on
    ylim([-30 0])
end

function [output_args,ERL,min_ERL]=TDR_ERL_Processing(output_args,OP,package_testcase_i,chdata,param)

%Fill TDR data
if package_testcase_i == 1
    if OP.TDR
        output_args.Z11est=chdata(1).TDR11.avgZport;
        if ~param.FLAG.S2P
            output_args.Z22est=chdata(1).TDR22.avgZport;
        else
            output_args.Z22est=[];
        end
        if OP.AUTO_TFX
            output_args.tfx_estimate=param.tfx(2);% 11/03/2021 RIM added for nbx estimate
        else
            output_args.tfx_estimate=[];
        end
    else
        output_args.Z11est=[];
        output_args.Z22est=[];
        output_args.tfx_estimate=[];
    end
end

% Process ERL
if package_testcase_i == 1
    if OP.ERL
        output_args.ERL11=chdata(1).TDR11.ERL;
        if ~param.FLAG.S2P
            output_args.ERL22=chdata(1).TDR22.ERL;
        else
            output_args.ERL22=[];
        end
        %                 output_args.ERL11RMS=chdata(1).TDR11.ERLRMS;
        %                  if ~param.FLAG.S2P,output_args.ERL22RMS=chdata(1).TDR22.ERLRMS;
    else
        output_args.ERL11=[];
        output_args.ERL22=[];
    end
end
if OP.ERL
    if OP.TDR_W_TXPKG
        min_ERL=output_args.ERL22;
        ERL= [ nan output_args.ERL22 ];
    else
        if ~isfield(output_args,'ERL22') || isempty(output_args.ERL22)
            min_ERL=output_args.ERL11;
            ERL= [ output_args.ERL11 nan ];
        else
            min_ERL=min(output_args.ERL11,output_args.ERL22);
            ERL= [ output_args.ERL11 output_args.ERL22 ];
        end
    end
    output_args.ERL=min_ERL;
else
    min_ERL=[];
    ERL= [];
    output_args.ERL=[];
end
if OP.ERL_ONLY
    if OP.BREAD_CRUMBS
        output_args.OP=OP;
        output_args.param=param;
        output_args.chdata=chdata;
        %This seems to be the intent of setting fom_result.ran=0.  Add it
        %to output_args so there is a fom_result field.
        fom_result.ran=0;
        output_args.fom_result=fom_result;
    end
    output_args.Z_t=param.Z_t;
    fileset_str=str2csv({chdata(1).base});
    output_args.file_names=sprintf('"%s"', fileset_str);
    if OP.DISPLAY_WINDOW
        savefigs(param, OP);
    end
end
function [impulse_response, p1_ctle, p2_ctle, z_ctle] = TD_CTLE(ir_in, fb, f_z, f_p1, f_p2, kacdc_dB, oversampling)
%% Equation 93A-22 implemented in z-space and applied to the impulse response.
p1_ctle = -2*pi*f_p1;
p2_ctle = -2*pi*f_p2;
z_ctle = -2*pi*f_z*10^(kacdc_dB/20);
k_ctle = -p2_ctle;
bilinear_fs = 2*fb*oversampling;
p2d = (1+p2_ctle/bilinear_fs)./(1-p2_ctle/bilinear_fs);
p1d = (1+p1_ctle/bilinear_fs)./(1-p1_ctle/bilinear_fs);
zd = (1+z_ctle/bilinear_fs)./(1-z_ctle/bilinear_fs);
% kd = (bilinear_fs-z_ctle)/((bilinear_fs-p1_ctle)*(bilinear_fs-p2_ctle));
% allow for different pole zeros RIM 9-29-2015
kd = (bilinear_fs-z_ctle)./((bilinear_fs-p1_ctle)*(bilinear_fs-p2_ctle))*f_p1/f_z;
B_filt =k_ctle*kd*poly([zd, -1]);
A_filt=poly([p1d, p2d]);
impulse_response=filter(B_filt,A_filt,ir_in);

function  [chdata, param, SDDch,SDDp2p] = TD_FD_fillin(param, OP, chdata)
Over_sample=2;
num_files=length(chdata);
for i=1:num_files
    V=chdata(i).uneq_pulse_response;
    T=chdata(i).t;
    dt=T(2)-T(1);
    f=0:1/max(T):1/dt;
    chdata(i).faxis=f;
    f75=find(f >= param.fb*.75,1,'first');
    fnq=find(f >= param.fb*.5,1,'first');
    chdata(i).fmaxi = length(f);
    chdata(i).faxis = f;
    UI=param.ui; % unit interval
    M=param.samples_per_ui; % sample per UI
    N_v=param.N_v; % number of UI for Vf determination
    
    % filters
    H_bt=Bessel_Thomson_Filter(param,f,OP.Bessel_Thomson);
    H_bw=Butterworth_Filter(param,f,OP.Butterworth);
    H_ftr=H_bw.*H_bt;
    H_ftr=H_ftr(:);
    % fd of PR
    prr=sin(f*UI*pi)./(f*UI*pi); %nremoving sinc function to avoid using sig proc toolbox
    prr = prr(:);
    if f(1)==0, prr(1)=1; end %remove NaN
    fd=fft(V);
    fd=fd(1:floor(length(fd)/2)); % un process freq domain response
    
    %% get Vf
    shifting_vector=kron(ones(1,200) ,[ 1  zeros(1,M-1) ]) ;
    step_response=filter(V,1, shifting_vector);
    Vf=step_response(end);
    STEP=timeseries(step_response(1:length(shifting_vector)),T(1:length(shifting_vector)));
    %%
    % ILest=20.*log10(abs(fd(1:f75)/Vf/M/Over_sample))  -   20.*log10(abs(prr(1:f75))) - 20*log10(abs(squeeze(H_ftr(1:f75)))) ;  %removing db function to avoid using sig proc toolbox
    % figure
    % plot(f(1:f75),ILest)
    
    IL_conv=fd(1:f75)/Vf/M/Over_sample  ./  prr(1:f75) ./ H_ftr(1:f75) ;  %removing db function to avoid using sig proc toolbox
    % set same variables as get_s4p_files
    IL_fields={'sdd12_raw' 'sdd21_raw' 'sdd12_orig' 'sdd21_orig' 'sdd12' 'sdd21' 'sdd21p' 'sdd21f'};
    Zero_fields={'sdd22_raw' 'sdd11_raw' 'sdc12_raw' 'sdc21_raw' 'sdc22_raw' 'sdc11_raw' ...
        'sdd11_orig' 'sdd22_orig' 'sdd11' 'sdd22' 'sdc12' 'sdc21' 'sdc11' 'sdc22' 'sdc21p'};
    zero_vector=zeros(length(IL_conv),1);
    for j=1:length(IL_fields)
        chdata(i).(IL_fields{j})=IL_conv;
    end
    for j=1:length(Zero_fields)
        chdata(i).(Zero_fields{j})=zero_vector;
    end
    
    if i==1
        SDDch(:,1,2)=chdata.sdd12_raw;
        SDDch(:,2,1)=chdata.sdd21_raw;
        SDDch(:,1,1)=chdata.sdd11_raw;
        SDDch(:,2,2)=chdata.sdd22_raw;
        SDDp2p= zeros(length(IL_conv),1);
    end
    chdata(i).TX_RL=[];
    chdata(i).TDR11=[];
    chdata(i).PDTR11=[];
    chdata(i).TDR22=[];
    chdata(i).PDTR22=[];  
end


function H_tw=Tukey_Window(f,param,fr,fb)
% RIM 05/26/2022 added optional fr and fb
% fr is the start of the raised cosine window
% fb is the end of the raised cosine window
if ~exist('fr','var') && ~exist('fb','var')
    fb=param.fb;
    fr=param.f_r*param.fb;
end
fperiod=2*(fb-fr);
H_tw = [ ones(1,length(f(f<fr))) ...
    0.5*cos(2*pi*(f( f>=fr & f<=fb)  -fb ) /fperiod-pi)+.5 ...
    zeros(1,length(f(f>fb)) )];
H_tw=H_tw(1:length(f));
if 0
    plot(f/1e9,H_tw)
end



%% moved output control to functions
function [H_TxFFE] = Tx_FFE_Filter(varargin)
% Author: Richard Mellitz
% Date: 7/29/2022
% generate FD tx ffe system function
% varagins...
% param - stucture 
%     param.fb baud rate
%     param.Tx_FFE Tx FFE coef
% f - freq array
% Use_Tx_FFE = flag to use or not
% H_TxFFE is system function for Tx_FFE
db = @(x) 20*log10(abs(x));
[param,varargin]=varargin_extractor(varargin{:});
[f,varargin]=varargin_extractor(varargin{:});
[Use_Tx_FFE,varargin]=varargin_extractor(varargin{:});
if isempty(Use_Tx_FFE)
    Use_Tx_FFE=0;
end
if isempty(param)
    param.fb=106.25e9;
    Tx_FFE=[1 ];
else
    if ~isfield(param, 'Pkg_TXFFE_preset')
        Tx_FFE=[  1 ];
    else
        Tx_FFE=param.Pkg_TXFFE_preset;
    end
end
if isempty(f)
    f=0:10e6:param.fb;
end


if Use_Tx_FFE ~=0
    [mcur,icur] = max(Tx_FFE);
    H_TxFFE=zeros(1,length(f));
    for ii=1:length(Tx_FFE)
        H_TxFFE = Tx_FFE(ii).*exp(-1j*2*pi*(ii-icur).*f/param.fb)+H_TxFFE;
    end
else
    H_TxFFE=ones(1,length(f));
end
% figure (1102320)
% plot(f/1e9,db(H_TxFFE))
% hold on
function C= WIENER_HOPF_MMSE(vsampled ,param, OP , chdata, txffe, Noise_XC,ivs) 
cmx=param.RxFFE_cmx;
cpx=param.RxFFE_cpx;
% do this early on so we can reuse the old code
% to be replaced with MMSE function from healey...
if param.N_bg ~=0 % must be floating taps
    cpx=param.N_bmax; % N_f in spreadsheet
end
num_taps=cmx+cpx+1;
ndfe=param.ndfe;
spui=param.samples_per_ui;
%% Start of WIENER-HOPF MMSE EQ code
                R_n = zeros(num_taps,num_taps);
                R_xt = zeros(num_taps,num_taps);

                if OP.Do_Colored_Noise
                    % Form Noise Autocorrelation matrix
                    Noise_XC = reshape (Noise_XC,1,[]);
                    len = length(Noise_XC);
                    if len < num_taps, Noise_XC = [Noise_XC,zeros(1,num_taps-len)]; end
                    Noise_XC = Noise_XC(1:num_taps);
                	R_n = toeplitz (Noise_XC,Noise_XC);
                end
                %% Calculate Cross Talk Correlation matrix at T intervals.
                if OP.Do_XT_Noise
                    % Calculate variance of Tx signal based on +/-1 outer limits
                    Tx_sigma = sqrt(  (param.levels^2-1)/(3*(param.levels-1)^2) );
                    for jj = 2:length(chdata)
                        if isequal(chdata(jj).type, 'FEXT') || isequal(chdata(jj).type, 'NEXT')
                            if isequal(chdata(jj).type, 'FEXT')
                                % len = length(chdata(jj).ctle_imp_response);
                                % ch_imp = ifft( fft(reshape(chdata(jj).ctle_imp_response,1,[])) .* fft(reshape(upsample(txffe, param.samples_per_ui),1,[]),len) );
                                ch_imp= FFE( txffe , cmx,param.samples_per_ui, chdata(jj).ctle_imp_response).';
                                ch_imp = filter(ones(param.samples_per_ui, 1), 1, ch_imp);
                            elseif isequal(chdata(jj).type, 'NEXT')
                                ch_imp = chdata(jj).ctle_imp_response;
                                ch_imp = filter(ones(param.samples_per_ui, 1), 1, ch_imp);
                            end
                            norms = zeros(1,spui);
                            for ii = 1:spui
                                norms(ii) = norm(ch_imp(ii:spui:end));
                            end
                            % Pick out sampling phase with largest noise contribution
                            [~,cursor] = max(norms);
                            sub_sample_ch = ch_imp(cursor:spui:end);
                            xc = xcorr(sub_sample_ch,num_taps)* Tx_sigma^2;
                            xc = xc(num_taps+1:end);
                            xc = xc(1:num_taps);
                            R = toeplitz (xc,xc);
                            R_xt = R_xt + R;
                        end
                    end
                end
                %% Noise + Cross Talk contribution to R matrix
                R_n_xc = zeros(num_taps+ndfe);
                R_n_xc(1:num_taps,1:num_taps) = R_n + R_xt ;
                %% For least means squares, we want to solve
                %
                %  ro =  |Ryy  Ryx| * w
                %        |Rxy  Rxx|
                % see Cioffi chapter 3, 3.7.3

                himp = vsampled;
                RefTap = cmx+1;
                Signal_Variance = mean ( [-1, -1/3, 1/3, 1].^2 );
                Ryy.r = [1:num_taps];
                Ryy.c = [1:num_taps];
                Ryx.r  = 1:num_taps;
                Ryx.c  = num_taps + (1:ndfe);
                Rxy.r  = num_taps + (1:ndfe);
                Rxy.c  = 1:num_taps;
                Rxx.r  = num_taps+(1:ndfe);
                Rxx.c  = num_taps+(1:ndfe);
                himp = reshape (himp,1,[]);
                %% ro is simply the channel response reversed in time
                himp_lr = fliplr(himp); % make sure himp has enough pre/post cursors, e.g. numtaps+ndfe
                himp_lr = [zeros(1,num_taps+ndfe), himp_lr, zeros(1,num_taps+ndfe)];
                [~,pk] = max(himp_lr);
                r_indx = (1:num_taps) - RefTap;
                ro = [himp_lr(pk+r_indx),  zeros(1,ndfe)].';
                ro = ro*Signal_Variance;
                %% Setup up the covariance matrix
                R = zeros(num_taps+ndfe);
                % Form Ryy
                % Note: important to use whole impulse response
                % not just the part that spans the FFE.
                [ryy,lags] = xcorr(himp,himp, num_taps-1);
                R(Ryy.r,Ryy.c) = toeplitz (ryy(num_taps:end));

                % Form Rxx
                R(Rxx.r,Rxx.c) = diag(ones(1,ndfe));

                % Form Ryx columns
                Ryx_indxs = (1:num_taps)-1;
                for jj = 0:ndfe-1
                    %             R(Ryx.r,Ryx.c(jj+1) ) = himp_lr(pk+1 + Ryx_indxs -jj ).';
                    R(Ryx.r,Ryx.c(jj+1) ) = himp_lr(pk-cmx-1 + Ryx_indxs -jj ).';
                end
                % Form Rxy rows
                R(Rxy.r,Rxy.c ) = R(Ryx.r,Ryx.c )';

                % add in Signal Variance
                R = R*Signal_Variance;
                Rtmp = R;
                % Add in Xt and colored noise terms
                R = R + R_n_xc;

                % SNR = 25 dB
                SNR = OP.FFE_SNR;
                Noise_var =   max(vsampled)^2*Signal_Variance/10^(SNR/10);
                R_noise = diag(ones(1,num_taps))*Noise_var;
                if ~OP.Do_Colored_Noise &&  ~OP.Do_XT_Noise
                    R(Ryy.r,Ryy.c) = R(Ryy.r,Ryy.c) + R_noise;
                end


                %% Solve for equalizer weights
                w = inv(R)*ro;
                C = w;
                %% Deal with 1st post Cursor DFE weight saturation
                %  ro = Rw by moving "saturated" weights over to the LHS
                DFE_h1_indx = num_taps+1;
                Indx_full = 1:length(C);
                ws = C;
                if ndfe>0 && abs(C(DFE_h1_indx)) > param.bmax(1)
                    rtmp = reshape (ro,[],1);
                    Rtmp = R;
              	    % Move saturated DFE weights over to left hand side of equation
                    ws = zeros (size(C));
                    ws (DFE_h1_indx) = sign(C(DFE_h1_indx))*param.bmax(1);
             	    rtmp = rtmp - Rtmp*ws;

             	    % and remove the corresponding column from R
              	    Rtmp(:,DFE_h1_indx) = [];
                    Indx_full (DFE_h1_indx) = [];
                    % now Rtmp isn't square so have to use the R'R trick
                    % Probably a little dicey "theoretically" because
                    % w = inv(R)*ro is already the mmse solution
                    % now we at doing a R'R operation, but hey
                    % seems to work
                    % Alternative, since R is now over specified, more rows than
                    % columns, one could try removing one of the DFE rows from the
                    % Rxy  Rxx portion of the R matrix.

              	    w_partial = inv(Rtmp'*Rtmp)*(Rtmp'*rtmp);
              	    ws (Indx_full,:) = w_partial;
                    C = ws;
                end
                % From Cioffi, Chapter 3
                var_ffe_dfe = Signal_Variance-ws.'*ro;
                SNR_Eq = 10*log10(Signal_Variance/var_ffe_dfe-1);

                %% Scale FFE gain to target output voltage
                Target_ouput = vsampled(ivs)*10^(param.current_ffegain/20);
                C = C*Target_ouput;
                C = C(1:num_taps);
            %% End MMSE dfe code
function Write_CSV(output_args,csv_file)

items = fieldnames(output_args);
item_value_strings = cell(size(items));
for field_id=1:length(items)
    field_name=items{field_id};
    field_value=output_args.(field_name);
    if isstruct(output_args.(field_name))
        field_value='struct';
    end
    if ischar(field_value)
        item_value_strings{field_id}=field_value;
    elseif isempty(field_value)
        item_value_strings{field_id}='';
    elseif numel(field_value)==1
        item_value_strings{field_id}=num2str(field_value);
    else
        item_value_strings{field_id}=sprintf('"%s"', mat2str(field_value));
    end
end

header_string = str2csv(items);
data_string = str2csv(item_value_strings);
fid = fopen(csv_file,'w');
fprintf(fid,'%s\n', header_string);
fprintf(fid,'%s\n', data_string);
fclose(fid);
function [ s11out, s12out, s21out, s22out ] = add_brd(chdata, param, OP)
%% Used in Clause 92 for adding board trace between TP0 and TP2

switch chdata.type
    case 'THRU'
        z_bp_tx = param.z_bp_tx;
        z_bp_rx = param.z_bp_rx;
    case 'NEXT'
        z_bp_tx = param.z_bp_rx;
        z_bp_rx = param.z_bp_next;
    case 'FEXT'
        z_bp_tx = param.z_bp_fext;
        z_bp_rx = param.z_bp_rx;
end
% Same cap on each tx and rx three is a data stratue for bifrucation but
% logic no implemented here RIM 06/28/2019
zref=param.Z0;
c1=param.C_0;
c2=param.C_1;
f=chdata.faxis;
f(f<eps)=eps;
s11pad1t= -1i*2*pi.*f*c1(1)*zref./(2+1i*2*pi.*f*c1(1)*zref);
s21pad1t= 2./(2+1i*2*pi.*f*c1(1)*zref);
s11pad2t= -1i*2*pi.*f*c2(1)*zref./(2+1i*2*pi.*f*c2(1)*zref);
s21pad2t= 2./(2+1i*2*pi.*f*c2(1)*zref);
[ s11tx, s12tx, s21tx, s22tx ] = synth_tline(chdata.faxis, param.brd_Z_c(1), param.Z0, param.brd_gamma0_a1_a2, param.brd_tau, z_bp_tx);
% add Tx caps
[s11tx, s12tx, s21tx, s22tx  ]= ...
    combines4p(  s11pad1t, s21pad1t, s21pad1t, s11pad1t, s11tx, s12tx, s21tx, s22tx   );
[s11tx, s12tx, s21tx, s22tx  ]= ...
    combines4p(   s11tx, s12tx, s21tx, s22tx,s11pad2t, s21pad2t, s21pad2t, s11pad2t  );


s11pad1r= -1i*2*pi.*f*c1(2)*zref./(2+1i*2*pi.*f*c1(2)*zref);
s21pad1r= 2./(2+1i*2*pi.*f*c1(2)*zref);
s11pad2r= -1i*2*pi.*f*c2(2)*zref./(2+1i*2*pi.*f*c2(2)*zref);
s21pad2r= 2./(2+1i*2*pi.*f*c2(2)*zref);
[ s11rx, s12rx, s21rx, s22rx ] = synth_tline(chdata.faxis, param.brd_Z_c(2), param.Z0, param.brd_gamma0_a1_a2, param.brd_tau, z_bp_rx);
% add Rx caps
[s11rx, s12rx, s21rx, s22rx  ]= ...
    combines4p(  s11pad2r, s21pad2r, s21pad2r, s11pad2r, s11rx, s12rx, s21rx, s22rx   );
[s11rx, s12rx, s21rx, s22rx  ]= ...
    combines4p(   s11rx, s12rx, s21rx, s22rx,s11pad1r, s21pad1r, s21pad1r, s11pad1r  );


switch OP.include_pcb
    case 1
        [ s11out1, s12out1, s21out1, s22out1 ]=combines4p(  s11tx, s12tx, s21tx, s22tx, chdata.sdd11_raw, chdata.sdd12_raw, chdata.sdd21_raw, chdata.sdd22_raw );
        [ s11out, s12out, s21out, s22out ]=combines4p(  s11out1, s12out1, s21out1, s22out1,  s11rx, s12rx, s21rx, s22rx);
    case 2
        [ s11out, s12out, s21out, s22out ]=combines4p( chdata.sdd11_raw, chdata.sdd12_raw, chdata.sdd21_raw, chdata.sdd22_raw, s11rx, s12rx, s21rx, s22rx);
end

function [ s11out, s12out, s21out, s22out ] = add_brdorig(chdata, param, OP)
% Used in Clause 92 for adding board trace between TP0 and TP2
switch chdata.type
    case 'THRU'
        z_bp_tx = param.z_bp_tx;
    case 'NEXT'
        z_bp_tx = param.z_bp_next;
    case 'FEXT'
        z_bp_tx = param.z_bp_fext;
end

[ s11tx, s12tx, s21tx, s22tx ] = synth_tline(chdata.faxis, param.brd_Z_c(1), param.Z0, param.brd_gamma0_a1_a2, param.brd_tau, z_bp_tx);
[ s11rx, s12rx, s21rx, s22rx ] = synth_tline(chdata.faxis, param.brd_Z_c(2), param.Z0, param.brd_gamma0_a1_a2, param.brd_tau, param.z_bp_rx);

switch OP.include_pcb
    case 1
        [ s11out1, s12out1, s21out1, s22out1 ]=combines4p(  s11tx, s12tx, s21tx, s22tx, chdata.sdd11_raw, chdata.sdd12_raw, chdata.sdd21_raw, chdata.sdd22_raw );
        [ s11out, s12out, s21out, s22out ]=combines4p(  s11out1, s12out1, s21out1, s22out1,  s11rx, s12rx, s21rx, s22rx);
    case 2
        [ s11out, s12out, s21out, s22out ]=combines4p( chdata.sdd11_raw, chdata.sdd12_raw, chdata.sdd21_raw, chdata.sdd22_raw, s11rx, s12rx, s21rx, s22rx);
end

function [hisi,tap_coef,hisi_ref]=applyDFEbk(hisi,hisi_ref,idx,tap_bk,curval,bmaxg,dfe_delta)
% hrem=applyDFEbk(hsis,idx,tap_bk,bmaxg)
% applying a bank of DFE at desired location
% hisi: waveform with cursor values;
% idx: starting index of the bank;
% tap_bk: number of taps in bank;
% bmaxg: maxmum coefficient in bank;

if nargin<6
    dfe_delta=0;
end

rng=idx:idx+tap_bk-1;
flt_curval=hisi(rng);

%floor(abs(floatingcursors/sbr(cursor_i))./param.dfe_delta).*param.dfe_delta.*sign(floatingcursors)*sbr(cursor_i);

if dfe_delta~=0
    flt_curval_q=floor(abs(flt_curval/curval)./dfe_delta) .* ...
        dfe_delta.*sign(flt_curval)*curval;
else
    flt_curval_q=hisi(rng);
end

tap_coef=min(abs(flt_curval_q/curval),bmaxg).*sign(flt_curval_q);
hisi(rng)= hisi(rng) - curval*tap_coef;
hisi_ref(rng)=0;

%AJG021820
function [ a] = bessel( n )
% bessel polynomial
for ii= 0:n
    a(ii+1) = factorial(2*n-ii) /  (2^(n-ii)*factorial(ii)*factorial(n-ii));
end


function [delay_sec, delay_idx]= calculate_delay_CausalityEnforcement(freq, sdd21, param, OP)
% History:
% 1. 14th October, 2021 (Intial release)

% Definition:
% This function captures the channel delay through the time domain using causality enforcement. 
% Following are the steps being followed.
% Step 1. Cascade negative frequencies
% Step 2. Extract magnitude
% Step 3. IFFT of the magnitude
% Step 4. Multiply by the sign(t)
% Step 5. Calculate the phase of the 1j*causal_phase
% Step 6. casual_function= |original|*exp(-1j*causal_phase)
% Step 7. f-domain to t-domain pulse response
% Step 8. Calculate the delay

% Author:
% Hansel Dsilva (dsilvahansel@gmail.com or hanseldsilva@achronix.com)

% Reference:
% 1] "IEEE Standard for Electrical Characterization of Printed Circuit Board and Related Interconnects at Frequencies up to 50 GHz," in IEEE Std 370-2020 , vol., no., pp.1-147, 8 Jan. 2021, doi: 10.1109/IEEESTD.2021.9316329.

% Input:
% freq          %frequency in hertz (odd number points)
% sdd21          %insertion loss in complex (odd number points)
% param          %COM native structure passed
% OP          %COM native structure passed

% Output:
% delay_sec          %channel delay in seconds
% delay_idx          %channel delay index

if iscolumn(sdd21)
    sdd21= sdd21.';
end
if iscolumn(freq)
    freq= freq.';
end

%---start. Step 1. Cascade negative frequencies
% sdd21_conj= zeros(1, length(freq)+ length(freq) - 1);
sdd21_conj= [sdd21, conj(sdd21(end:-1:2))];%For some reason this only works 
% sdd21_conj = [real(sdd21(1)), sdd21(2:end-1), real(sdd21(end)), flipud(conj(sdd21(2:end-1)))];
%---end. Step 1. Cascade negative frequencies

%---start. Step 2. Extract magnitude
sdd21_mag_conj = real(log(abs(sdd21_conj)));
%---end. Step 2. Extract magnitude

%---start. Step 3. IFFT of the magnitude
sdd21_mag_time = ifft(sdd21_mag_conj);
%---end. Step 3. IFFT of the magnitude

%---start. Step 4. Multiply by the sign(t)
sdd21_mag_time(1:length(freq))= +1j*sdd21_mag_time(1:length(freq));
sdd21_mag_time(length(freq)+1:end)= -1j*sdd21_mag_time(length(freq)+1:end);
%---end. Step 4. Multiply by the sign(t)

%---start. Step 5. Calculate the phase of the original_signal and the 1j*causal_phase
sdd21_phase_causality_enforced = real(fft(sdd21_mag_time));
%---end. Step 5. Calculate the phase of the original_signal and the 1j*causal_phase

%---start. Step 6. casual_function= |original|*exp(-1j*causal_phase)
sdd21_causality_enforced= abs(sdd21_conj).*exp(-1j*sdd21_phase_causality_enforced);
sdd21_causality_enforced= sdd21_causality_enforced(1:length(freq));
%---end. Step 6. casual_function= |original|*exp(-1j*causal_phase)

%---start. Step 7. f-domain to t-domain pulse response
%------Note. Do not use s21_to_impulse() for we do not want to truncate
%--------- Extrapolation has been already done by the COM tool
freq_array= freq;
time_step= param.sample_dt;
fmax=1/time_step/2;
freq_step=(freq_array(3)-freq_array(2))/1;
fout=0:1/round(fmax/freq_step)*fmax:fmax;

ILin=sdd21;
IL=interp_Sparam(ILin,freq_array,fout, ...
    OP.interp_sparam_mag, OP.interp_sparam_phase,OP);
IL_nan = find(isnan(IL));
for in=IL_nan
    IL(in)=IL(in-1);
end
IL = IL(:);
% add padding for time steps
IL_symmetric = [real(IL(1)); IL(2:end-1); real(IL(end)); flipud(conj(IL(2:end-1)))];
sdd21_PR = filter(ones(1, param.samples_per_ui), 1, real(ifft(IL_symmetric)));%real(ifft(IL_symmetric));
clear IL IL_nan IL_symmetric

ILin=sdd21_causality_enforced;
IL=interp_Sparam(ILin,freq_array,fout, ...
    OP.interp_sparam_mag, OP.interp_sparam_phase,OP);
IL_nan = find(isnan(IL));
for in=IL_nan
    IL(in)=IL(in-1);
end
IL = IL(:);
% add padding for time steps
% IL_symmetric = [IL(1:end-1); 0; flipud(conj(IL(2:end-1)))];
% IL_symmetric = [IL(1:end); flipud(conj(IL(2:end-1)))];
IL_symmetric = [real(IL(1)); IL(2:end-1); real(IL(end)); flipud(conj(IL(2:end-1)))];
sdd21_causality_enforced_PR = filter(ones(1, param.samples_per_ui), 1, real(ifft(IL_symmetric)));%real(ifft(IL_symmetric));
clear IL IL_nan IL_symmetric

clear time_step fmax freq_step freq_array

freq_step=(fout(3)-fout(2))/1;
L= length(sdd21_PR);
t_base = (0:L-1)/(freq_step*L);
clear fout freq_step L
%---end. Step 7. f-domain to t-domain pulse response

%---start. Step 8. Calculate the delay
%------start. Remove the last 5% of the waveform for noise due to IFFT
sdd21_causality_enforced_PR_reduced= sdd21_PR(1:floor(end*95/100));
sdd21_PR_reduced= sdd21_causality_enforced_PR(1:floor(end*95/100));
%------end. Remove the last 5% of the waveform for noise due to IFFT

%---start. calculate the difference in index between the peaks
[~, peak_x_idx] = max(sdd21_causality_enforced_PR_reduced);
[~, peak_y_idx] = max(sdd21_PR_reduced);
peak_idx_difference = peak_x_idx - peak_y_idx;
%---end. calculate the difference in index between the peaks

if peak_idx_difference~=0
    search_bounds = min(peak_x_idx, peak_y_idx);%<----one may fix it to 1e3; using minimum of the peaks in assuming the other signal has zero delay
    error_value = length(sdd21_causality_enforced_PR_reduced);
    error_idx = 0;
    %     i= 1;
    for shift_value= peak_idx_difference-search_bounds:peak_idx_difference+search_bounds
        sdd21_PR_reduced_shifted = circshift(sdd21_PR_reduced,shift_value);
        current_error = norm(sdd21_PR_reduced_shifted-sdd21_PR_reduced);
        if (error_value > current_error)
            error_idx = shift_value;
            error_value = current_error;
        end
        %         error_idx_H(i)= error_idx;
        %         i= i+ 1;
    end
    %plot(error_idx_H);
    
    if error_idx==0
        error('An odd case has occured when calcualting the channel delay. Please contact the tool developer.');
    end
    
    delay_idx = error_idx;
else
    delay_idx= 1;
end

delay_sec= t_base(abs(delay_idx));


function [return_struct]= capture_RIL_RILN(chdata)
% History:
% 1. 12th April, 2019 (Intial release)
%
% 2. 11th October, 2021
%   - Details:
%       1] Revised the equation of RIL for missing conj(rho_port1) and conj(rho_port2).
%       2] Revised the selection criteria for the solution of the quadratic
%       equation in finding the reflection coefficient (rho).
%   - Impact:
%       => Zero impact in |RIL|, while impact on angle(RIL).
%   - Previous:
%       %---start. For passive networks the reflection coefficient should be less than one
%       if all(abs(solution_1(idx_start:end))< 1) && ~all(abs(solution_2(idx_start:end))< 1)
%           rho_port2= solution_1;
%        elseif all(abs(solution_2(idx_start:end))< 1) && ~all(abs(solution_1(idx_start:end))< 1)
%          rho_port2= solution_2;
%        else
%            rho_port2= solution_1;
%           %     error('Please contact the tool developer. It appears a odd case has appeared.');
%       end
%       %---end. For passive networks the reflection coefficient should be less than one
%
%       RIL= conj((1- rho_port1).*(sqrt(abs(1- rho_port1.*conj(rho_port1)))./abs(1-rho_port1)  )).*(  Sdd21.*(1- rho_port2.*Sdd22)+ (Sdd22- conj(rho_port2)).*rho_port2.*Sdd21  )./ ( ((1- rho_port2).*(sqrt(abs(1- rho_port2.*conj(rho_port2)))./abs(1-rho_port2)  )) .* (  -rho_port1.*rho_port2.*Sdd21.*Sdd12+  (1- rho_port2.*Sdd22).*(1- rho_port1.*Sdd11)  ) );
%   - Change:
%       %---start. Given the real part of the impedance is to be positive
%       Z_solution_1= (solution_1*conj(SCH.Impedance)+ SCH.Impedance)./(1-solution_1);
%       Z_solution_2= (solution_2*conj(SCH.Impedance)+ SCH.Impedance)./(1-solution_2);
% 
%       rho_port2= zeros(length(solution_1), 1);
%        for solution_idx= 1:length(solution_1)
%          if real(Z_solution_1(solution_idx))>0 && real(Z_solution_2(solution_idx))<=0
%              rho_port2(solution_idx, 1)= solution_1(solution_idx);
%           elseif real(Z_solution_2(solution_idx))>0 && real(Z_solution_1(solution_idx))<=0
%              rho_port2(solution_idx, 1)= solution_2(solution_idx);
%          else
%               error('An odd case has occured. Please contact the tool developer.');        
%           end
%       end
%       %---end. Given the real part of the impedance is to be positive
%       RIL= conj((1- conj(rho_port1)).*(sqrt(abs(1- rho_port1.*conj(rho_port1)))./abs(1-rho_port1)  )).*(  Sdd21.*(1- rho_port2.*Sdd22)+ (Sdd22- conj(rho_port2)).*rho_port2.*Sdd21  )./ ( ((1- conj(rho_port2)).*(sqrt(abs(1- rho_port2.*conj(rho_port2)))./abs(1-rho_port2)  )) .* (  -rho_port1.*rho_port2.*Sdd21.*Sdd12+  (1- rho_port2.*Sdd22).*(1- rho_port1.*Sdd11)  ) );

% Definition:
% This function captures the reflectionless insertion loss (RIL) and reflective insertion loss nois (RILN) for any arbitary S-parameter

% Author:
% Hansel Dsilva (dsilvahansel@gmail.com or hanseldsilva@achronix.com)
% Acknowledgement: Adam Gregory, Richard Mellitz, Beomtaek Lee and Amit Kumar.

% This function has been shared by Hansel for others to evaluate for reflectionless insertion loss (RIL) and reflective insertion loss nois (RILN) for any arbitary S-parameter.

% Reference:
% 1] H. Dsilva et al., "Finding Reflective Insertion Loss Noise and Reflectionless Insertion Loss," 2020 DesignCon.
% 2] H. Dsilva, A. Jain, J. Sasikala and A. Kumar, "Novel Signal Integrity Application of Power Wave Scattering Matrix theory," 2019 IEEE MTT-S International Microwave and RF Conference (IMARC), 2019, pp. 1-7, doi: 10.1109/IMaRC45935.2019.9118701.

% Input:
% 1] SCH: S-matrix structure
% SCH.Frequencies= faxis;
% SCH.Parameters(1,1,:)= sdd11;
% SCH.Parameters(2,2,:)= sdd22;
% SCH.Parameters(1,2,:)= sdd12;
% SCH.Parameters(2,1,:)= sdd21;
% SCH.NumPorts= 2;
% SCH.Impedance= 100;

% Output: Struct returned has the following,
% return_struct.RIL          %Reflectionless Insertion Loss as a complex number
% return_struct.RIL_dB          %Reflectionless Insertion Loss in decibel
% return_struct.RILN          %Reflective Insertion Loss Noise as a complex number
% return_struct.RILN_dB          %Reflective Insertion Loss Noise in decibel
% return_struct.Z_port1          %Frequency dependent, complex values of impedance for termination of port 1
% return_struct.Z_port2          %Frequency dependent, complex values of impedance for termination of port 2
% return_struct.rho_port1          %Frequency dependent, complex values of the reflection coefficient of port 1
% return_struct.rho_port2          %Frequency dependent, complex values of reflection coefficient of port 2
% return_struct.freq          %Frequency axis

SCH.Parameters(1,1,:)=chdata(1).sdd11_orig;
SCH.Parameters(2,2,:)=chdata(1).sdd22_orig;
SCH.Parameters(1,2,:)=chdata(1).sdd12_orig;
SCH.Parameters(2,1,:)=chdata(1).sdd21_orig;
SCH.Frequencies=chdata(1).faxis;
SCH.Impedance=100;%<------------------in a future release, may want to parameterize this based on the read .sNp
SCH.NumPorts= 2;

%---start. allowed is only a 2-port network having a transmitter and receiver
if size(SCH.Parameters, 1)~=2
	fprintf('Size of given S matrix is %d.', size(SCH.Parameters, 3) );
	error('Allowed is only a 2-port network having a transmitter and receiver.');
end
%---end. allowed is only a 2-port network having a transmitter and receiver

%---start. do not include the DC point given sinusoidals at DC are not
%defined
if SCH.Frequencies(1)==0
    idx_start= 2;
else
    idx_start= 1;
end
%---end. do not include the DC point given sinusoidals at DC are not
%defined

Sdd11= squeeze(SCH.Parameters(1,1,idx_start:end));
Sdd12= squeeze(SCH.Parameters(1,2,idx_start:end));
Sdd21= squeeze(SCH.Parameters(2,1,idx_start:end));
Sdd22= squeeze(SCH.Parameters(2,2,idx_start:end));

a= -Sdd22+ Sdd11.*Sdd22.*conj(Sdd11)- Sdd21.*Sdd12.*conj(Sdd11);
b= 1+ Sdd22.*conj(Sdd22)+...
    Sdd11.*Sdd22.*conj(Sdd21).*conj(Sdd12)-...
    Sdd21.*Sdd12.*conj(Sdd21).*conj(Sdd12)-...
    Sdd11.*Sdd22.*conj(Sdd11).*conj(Sdd22)+...
    Sdd12.*Sdd21.*conj(Sdd11).*conj(Sdd22)-...
    Sdd11.*conj(Sdd11);
c= -conj(Sdd22)-...
    Sdd11.*conj(Sdd21).*conj(Sdd12)+...
    Sdd11.*conj(Sdd11).*conj(Sdd22);

solution_1= (-b+sqrt((b.^2)-(4*a.*c)))./(2*a);
solution_2= (-b-sqrt((b.^2)-(4*a.*c)))./(2*a);

clear a b c

%---start. Given the real part of the impedance is to be positive
Z_solution_1= (solution_1*conj(SCH.Impedance)+ SCH.Impedance)./(1-solution_1);
Z_solution_2= (solution_2*conj(SCH.Impedance)+ SCH.Impedance)./(1-solution_2);

rho_port2= zeros(length(solution_1), 1);
for solution_idx= 1:length(solution_1)
    if real(Z_solution_1(solution_idx))>0 && real(Z_solution_2(solution_idx))<=0
        rho_port2(solution_idx, 1)= solution_1(solution_idx);
    elseif real(Z_solution_2(solution_idx))>0 && real(Z_solution_1(solution_idx))<=0
        rho_port2(solution_idx, 1)= solution_2(solution_idx);
    else
        error('An odd case has occured. Please contact the tool developer.');        
    end
end
%---end. Given the real part of the impedance is to be positive

rho_port1= conj(Sdd11+ ((rho_port2.*Sdd21.*Sdd12)./(1-(rho_port2.*Sdd22))));

%---start. calculate the equivalent port impedance
Z_port1= (rho_port1*conj(SCH.Impedance)+ SCH.Impedance)./(1-rho_port1);
Z_port2= (rho_port2*conj(SCH.Impedance)+ SCH.Impedance)./(1-rho_port2);
%---end. calculate the equivalent port impedance


% %---start. The reflectionless insertion loss is the insertion loss corresponding
% %to zero reflections.
RIL= conj((1- conj(rho_port1)).*(sqrt(abs(1- rho_port1.*conj(rho_port1)))./abs(1-rho_port1)  )).*(  Sdd21.*(1- rho_port2.*Sdd22)+ (Sdd22- conj(rho_port2)).*rho_port2.*Sdd21  )./ ( ((1- conj(rho_port2)).*(sqrt(abs(1- rho_port2.*conj(rho_port2)))./abs(1-rho_port2)  )) .* (  -rho_port1.*rho_port2.*Sdd21.*Sdd12+  (1- rho_port2.*Sdd22).*(1- rho_port1.*Sdd11)  ) );
%---end. The reflectionless insertion loss is the insertion loss corresponding
%to zero reflections.

%---start. Calculate RILN (Reflective Insertion Loss Noise)
RILN= RIL- Sdd21;
RILN_dB= 20*log10(abs(Sdd21))- 20*log10(abs(RIL));
%---end. Calculate RILN (Reflective Insertion Loss Noise)

%---start.  preparing returns struct
return_struct.RIL= RIL;		%Reflectionless Insertion Loss as a complex number
return_struct.RIL_dB= 20*log10(abs(RIL));		%Reflectionless Insertion Loss in decibel
return_struct.RILN= RILN;		%Reflective Insertion Loss Noise as a complex number
return_struct.RILN_dB= RILN_dB;		%Reflective Insertion Loss Noise in decibel
return_struct.Z_port1= Z_port1;		%Frequency dependent, complex values of impedance for termination of port 1
return_struct.Z_port2= Z_port2;		%Frequency dependent, complex values of impedance for termination of port 2
return_struct.rho_port1= rho_port1;		%Frequency dependent, complex values of the reflection coefficient of port 1
return_struct.rho_port2= rho_port2;		%Frequency dependent, complex values of reflection coefficient of port 2
return_struct.freq= SCH.Frequencies(idx_start:end);		%Frequency axis
%---end. preparing returns struct
function [noise_bottom,noise_top]=cdf_to_ber_contour(cdf,specBER)


%For the given BER, find the top & bottom voltage level in the CDF

%for the top, just find the first index at the spec BER
nidx=find(cdf.y>specBER, 1, 'first');
noise_bottom = cdf.x(nidx);
%for top, flip the cdf.  need to operate on row vector for fliplr:  cdf.y(:)'
nidx=find(fliplr(cdf.y(:)')>specBER, 1, 'first');
%the true index without flipping
nidx=length(cdf.y)-nidx+1;
noise_top = cdf.x(nidx);
function p=comb_fct(p1, p2)
if p1.BinSize ~= p2.BinSize
    error('bin size must be equal')
end

p=p1;
p.BinSize=p1.BinSize;
%p.Min=p1.Min+p2.Min;
p.Min=min(p1.Min,p2.Min);
difsz=abs(p1.Min-p2.Min);
lp1=length(p1.y);
lp2=length(p2.y);
if p1.Min == p.Min
    p2.y(difsz+1:lp2+difsz)=p2.y;
    p2.y(1:difsz)=0;
    p2.y(lp2+difsz+1:lp1)=0;
elseif p2.Min == p.Min
    p1.y(difsz+1:lp1+difsz)=p1.y;
    p1.y(1:difsz)=0;
    p1.y(lp1+difsz+1:lp2)=0;
end
    % p.y=conv(p1.y, p2.y);
p.y=(p1.y+p2.y);
% p.y=p.y/sum(p.y);
% p.x =p.Min*p.BinSize:p.BinSize:-p.Min*p.BinSize;
p.x =(p.Min:-p.Min)*p.BinSize;	% modified by Yasuo Hidaka, 9/4/2016


function out_pdf=combine_pdf_same_voltage_axis(pdf1,pdf2)

if pdf1.BinSize ~= pdf2.BinSize
    error('bin size must be equal')
end

x1=pdf1.x;
y1=pdf1.y;
x2=pdf2.x;
y2=pdf2.y;
%find the pdf with a larger min, force it to have the same min, and insert
%probability = 0
min1=pdf1.x(1);
min2=pdf2.x(1);
shift_amount=round(abs(min1-min2)/pdf1.BinSize);
if min1<min2
    x2=[pdf1.x(1:shift_amount) pdf2.x];
    y2=[zeros(1,shift_amount) pdf2.y];
else
    x1=[pdf2.x(1:shift_amount) pdf1.x];
    y1=[zeros(1,shift_amount) pdf1.y];
end
%find the pdf with smaller max, force it to have the same max, and insert
%probability=0
L1=length(x1);
L2=length(x2);
Ldiff=abs(L1-L2);
if L1>L2
    out_x=x1;
    y2=[y2 zeros(1,Ldiff)];
else
    out_x=x2;
    y1=[y1 zeros(1,Ldiff)];
end
%now the 2 pdfs have the same voltage axis, add probabilities together
%renormalization is not handled here, so the output pdf will not have sum=1
%It is the responsibility of the calling function to handle renormalization
%if needed
out_y=y1+y2;
out_pdf.x=out_x;
out_pdf.y=out_y;
out_pdf.BinSize=pdf1.BinSize;
function [ s11out, s12out, s21out, s22out ] = combines4p( s11in1, s12in1, s21in1, s22in1, s11in2, s12in2, s21in2, s22in2)

%original method:
% s1=zeros(2,2,length(s11in1)); s2=s1; t3=s1;
% for i=1:length(s11in1)
%     s1(:,:,i)=[s11in1(i) s12in1(i); s21in1(i) s22in1(i) ];
%     s2(:,:,i)=[s11in2(i) s12in2(i); s21in2(i) s22in2(i) ];
% end
% t1=stot(s1);
% t2=stot(s2);
% for i=1:length(s11in1)
%     t3(:,:,i)=t1(:,:,i)*t2(:,:,i);
% end
% s3=ttos(t3);
% s11out=s3(1,1,:);
% s11out=transpose(s11out(:));
% s12out=s3(1,2,:);
% s12out=transpose(s12out(:));
% s21out=s3(2,1,:);
% s21out=transpose(s21out(:));
% s22out=s3(2,2,:);
% s22out=transpose(s22out(:));


%vectorized method:
s1(1,1,:)=s11in1;
s1(1,2,:)=s12in1;
s1(2,1,:)=s21in1;
s1(2,2,:)=s22in1;
s2(1,1,:)=s11in2;
s2(1,2,:)=s12in2;
s2(2,1,:)=s21in2;
s2(2,2,:)=s22in2;


N = (1-s1(2,2,:).*s2(1,1,:)) ;
s11out = s1(1,1,:)+(s1(1,2,:).*s1(2,1,:).*s2(1,1,:))./N ;
s12out = s1(1,2,:).*s2(1,2,:)./N ;
s21out = s2(2,1,:).*s1(2,1,:)./N;
s22out = s2(2,2,:)+(s2(1,2,:).*s2(2,1,:).*s1(2,2,:))./N ;

s11out=transpose(squeeze(s11out));
s12out=transpose(squeeze(s12out));
s21out=transpose(squeeze(s21out));
s22out=transpose(squeeze(s22out));
function p=conv_fct(p1, p2)
if p1.BinSize ~= p2.BinSize
    error('bin size must be equal')
end

p=p1;
%p.BinSize=p1.BinSize;
%p.Min=p1.Min+p2.Min;
p.Min=round(p1.Min+p2.Min);	% modified by Yasuo Hidaka, 9/4/2016
p.y=conv2(p1.y, p2.y);
%p.x =p.Min*p.BinSize:p.BinSize:-p.Min*p.BinSize;
%p.x =(p.Min:-p.Min)*p.BinSize;	% modified by Yasuo Hidaka, 9/4/2016
pMax=p.Min+length(p.y)-1;
p.x =(p.Min*p.BinSize:p.BinSize:pMax*p.BinSize);




function p=conv_fct_MeanNotZero(p1, p2)

if p1.BinSize ~= p2.BinSize
error('bin size must be equal')
end

p=p1;
%p.BinSize=p1.BinSize;
%p.Min=p1.Min+p2.Min;
p.Min=round(p1.Min+p2.Min);	% modified by Yasuo Hidaka, 9/4/2016
p.y=conv2(p1.y, p2.y);

%This is equivalent to (p.Min:p.Min+length(p.y)-1)*p.BinSize
%But it is faster to pre-multiply BinSize instead of multiplying the entire
%vector by BinSize
pMax=p.Min+length(p.y)-1;
p.x =(p.Min*p.BinSize:p.BinSize:pMax*p.BinSize);
function [cursor_i,no_zero_crossing,sbr_peak_i,zxi]=cursor_sample_index(sbr,param,OP,peak_search_range)

%IN:
%sbr = pulse response
%param = COM "param" struct
%OP = COM "OP" struct
%peak_search_range= a limited range to search for the peak (for speed up)
%       it is usually +/- 20 UI
%
%OUT:
%cursor_i = sampling location
%no_zero_crossing = flag that reveals if sampling is not possible.
%       When this function is called in optimize_fom, it signals to quit the current case
%sbr_peak_i = index of the pulse peak.  Note:  tens of thousands of calls to max(sbr) are very
%       time consuming, so saving the peak in one spot is advantageous
%zxi = zero crossing index (returned because RXFFE uses it)

no_zero_crossing=0;
%need to set cursor_i to empty in case no_zero_crossing flag is set
cursor_i=[];

%get peak of pulse and peak index
[max_of_sbr, sbr_peak_tmp]=max(sbr(peak_search_range));
sbr_peak_i=sbr_peak_tmp+peak_search_range(1)-1;


% initial guess at cursor location (t_s)  - based on approximate zero crossing
%limit search space for speed up
search_start=sbr_peak_i-4*param.samples_per_ui;
if search_start<1
    search_start=1;
end
%Find zero crossings
zxi = find(diff(sign(sbr(search_start:sbr_peak_i)-.01*max_of_sbr))>=1)+search_start-1;

%Note:  the original implementation of zxi:
% zxi = find(diff(sign(sbr-.01*max(sbr)))>=1);
% zxi = zxi(zxi<sbr_peak_i);
% zxi = zxi(sbr_peak_i - zxi < 4*param.samples_per_ui);
% The changes to limit search space and remember max(sbr) give 10x speed up
% A test case of 25k runs, reduced from 1.2s to 0.1s


if isempty(zxi)
    %if no zero crossing, the calling program must respond (since sample point will be empty)
    no_zero_crossing=1;
    return;
elseif length(zxi)>1
    %only need the last zero crossing
    zxi=zxi(end);
end
if param.ndfe==0
    max_dfe1=0;
else
    max_dfe1=param.bmax(1);
end
% adjust cursor_i to Solve equation 93A-25 %%
% Muller-Mueller criterion with DFE
mm_range = zxi+(0:2*param.samples_per_ui);
switch OP.CDR
    case 'Mod-MM'
        mm_metric = ...
            abs(sbr(mm_range+param.samples_per_ui)-max_dfe1*sbr(mm_range));
    otherwise % MM
        %MM is generally: first precursor = 0
        %but the actual requirement is for first precursor = first postcursor (after DFE is applied)
        %if first postcursor doesn't exceed max_dfe, then this equates to first precursor = 0
        %in cases where first postcursor exceeds max_dfe, the mismatch is balanced out so that
        %first precursor = first postcursor - max_dfe
        mm_metric = ...
            abs(sbr(mm_range-param.samples_per_ui) - max(sbr(mm_range+param.samples_per_ui)-max_dfe1*sbr(mm_range), 0));
end
[~, mm_cursor_offset] = min(mm_metric);

%cursor_i = the sample location
cursor_i = zxi+mm_cursor_offset-1;
function pdf=d_cpdf( binsize, values, probs)
%  p=cpdf(type, ...)
%
% CPDF is a probability mass function for discrete distributions or an
% approxmation of a PDF for continuous distributions.
%
% cpdf is internally normalized so that the sum of probabilities is 1
% (regardless of bin size).

% Internal fields:
% Min: *bin number* of minimum value.
% BinSize: size of PDF bins. Bin center is the representative value.
% Vec: vector of probabilities per bin.

%speed up for initializing empty pdf
if all(values==0)
    pdf.BinSize=binsize;
    pdf.Min=0;
    pdf.y=1;
    pdf.x=0;
    return;
end

if ~issorted(values)
    [values,si]=sort(values);
    probs=probs(si);
end
values=binsize*round(values/binsize);
t=(values(1):binsize:values(end));
pdf.Min=values(1)/binsize;
pdf.y=zeros(size(t));
for k=1:length(values)
    if k==1
        bin=1;
    elseif k==length(values)
        bin=length(t);
    else
        [UNUSED_OUTPUT, bin]=min(abs(t-values(k))); %#ok<ASGLU>
    end
    pdf.y(bin) = pdf.y(bin)+probs(k);
end

pdf.BinSize=binsize;
pdf.y=pdf.y/sum(pdf.y);
if any(~isreal(pdf.y)) || any(pdf.y<0)
    error('PDF must be real and nonnegative');
end
support=find(pdf.y);
pdf.y=pdf.y(support(1):support(end));
pdf.Min=pdf.Min+(support(1)-1);
pdf.x=(pdf.Min:-pdf.Min)*binsize;
function clip_output=dfe_clipper(input,max_threshold,min_threshold)

if isrow(input)
    max_threshold=max_threshold(:).';
    min_threshold=min_threshold(:).';
else
    max_threshold=max_threshold(:);
    min_threshold=min_threshold(:);
end

clip_output=input;
clip_output(input>max_threshold)=max_threshold(input>max_threshold);
clip_output(input<min_threshold)=min_threshold(input<min_threshold);


function [msg] = end_display_control(msg,param,OP,output_args,COM,min_ERL,ERL,VEO_mV,VEC_dB,threshold_DER,DISPLAY_WINDOW)
[ncases, mele]=size(param.z_p_tx_cases);
if mele ==2
    param.flex=2;
elseif mele==4
    param.flex=4;
elseif mele==1
    param.flex=1;
else
    error(springf('config file syntax error'))
end


if DISPLAY_WINDOW && ~OP.RX_CALIBRATION
    % display bathtub curves in one axis per test case.
    %             h=findall(0, 'Name', 'COM results');
    if ~exist('h','var')
        msgtext = cell(1, length(OP.pkg_len_select));
        msgcolor = 'g';
    else
        msgtext=get(findobj(h, 'type', 'text'), 'string');
        msgcolor = get(h, 'color');
        close(h); % will be recreated
    end
    msgctr=size(msgtext,1)+1;
    if ~OP.ERL_ONLY
        switch OP.PHY
            case 'C2M'
                if VEO_mV >= param.Min_VEO && VEO_mV <= param.Max_VEO
                    msg=sprintf('%s: EH = %.3f mV (pass)\n', ...
                        msg, VEO_mV);
                else
                    msg=sprintf('%s: EH = %.3f mV (FAIL)\n', ...
                        msg, VEO_mV);
                    msgcolor = 'r';
                end
                
                if VEC_dB <= param.VEC_pass_threshold
                    msg=sprintf('%s: VEC = %.3f dB (pass)\n', ...
                        (msg), VEC_dB);
                else
                    msg=sprintf('%s: VEC = %.3f dB (FAIL)\n', ...
                        (msg), VEC_dB);
                    msgcolor = 'r';
                end
            case 'C2C'
                if COM >= param.pass_threshold
                    %                         msgtext{package_testcase_i}=sprintf('%s: COM = %.3f dB (pass)\n', ...
                    %                             msg, COM);
                    msg=sprintf('%s: COM = %.3f dB (pass)\n', ...
                        msg, COM);
                else
                    %                         msgtext{package_testcase_i}=sprintf('%s: COM = %.3f dB (FAIL)\n', ...
                    %                             msg, COM);
                    msg=sprintf('%s: COM = %.3f dB (FAIL)\n', ...
                        msg, COM);
                    msgcolor = 'r';
                end
                % begin yasuo patch 3/18/2019
                msg=sprintf('%s: DER = %.3e at COM threshold \n', msg, threshold_DER);
                % end yasuo patch
            case 'C2Mcom'
                if VEO_mV >= param.Min_VEO && VEO_mV <= param.Max_VEO
                    msg=sprintf('%s: EH = %.3f mV (pass)\n', ...
                        msg, VEO_mV);
                else
                    msg=sprintf('%s: EH = %.3f mV (FAIL)\n', ...
                        msg, VEO_mV);
                    msgcolor = 'r';
                end
                
                if VEC_dB <= param.VEC_pass_threshold
                    msg=sprintf('%s: VEC = %.3f dB (pass)\n', ...
                        (msg), VEC_dB);
                else
                    msg=sprintf('%s: VEC = %.3f dB (FAIL)\n', ...
                        (msg), VEC_dB);
                    msgcolor = 'r';
                end
                if COM >= param.pass_threshold
                    %                         msgtext{package_testcase_i}=sprintf('%s: COM = %.3f dB (pass)\n', ...
                    %                             msg, COM);
                    msg=sprintf('%s: COM = %.3f dB (pass)\n', ...
                        msg, COM);
                else
                    %                         msgtext{package_testcase_i}=sprintf('%s: COM = %.3f dB (FAIL)\n', ...
                    %                             msg, COM);
                    msg=sprintf('%s: COM = %.3f dB (FAIL)\n', ...
                        msg, COM);
                    msgcolor = 'r';
                end
                % begin yasuo patch 3/18/2019
                msg=sprintf('%s: DER = %.3e at COM threshold \n', msg, threshold_DER);
                % end yasuo patch
        end
    end
    if OP.ERL
        if ~isempty(ERL)
            if min_ERL >= param.ERL_pass_threshold
                %                             msgtext{package_testcase_i+1}=[sprintf(' PASS ... ERL = %.3f dB\n', min_ERL)];
                msg=[sprintf('%s: PASS ... ERL = %.3f dB (%.3f dB,%.3f dB) \n',msg, min_ERL, ERL)];
            else
                %                             msgtext{package_testcase_i+1}=[ sprintf(' FAIL ... ERL = %.3f dB\n', min_ERL) ];
                msg=[ sprintf('%s: FAIL ... ERL = %.3f dB (%.3f dB,%.3f dB) \n',msg, min_ERL, ERL) ];
                msgcolor = 'r';
            end
        end
    end
    h=msgbox(msg, ['COM r' output_args.code_revision ' results']);
    set(h, 'color', msgcolor, 'tag', 'COM');
    movegui(h, 'center');
else % no windows
    %     display(['max noise at BER = ' num2str(peak_interference_at_BER)])
    %     display(['signal after eq = ' num2str(A_s/(param.levels-1))])
    if ~OP.ERL_ONLY
        switch OP.PHY
            case 'C2C'
                if COM >= param.pass_threshold
                    fprintf('%s <strong> PASS ... COM = %.3f dB</strong>\n', msg, COM);
                else
                    fprintf(2,'%s <strong> FAIL ... COM = %.3f dB</strong>\n',  msg, COM);
                end
                % begin yasuo patch 3/18/2019
                fprintf('%s DER = %.3e at COM threshold \n', msg, threshold_DER);
                % end yasuo patch
            case 'C2Mcom'
                if VEC_dB <= param.VEC_pass_threshold
                    fprintf('%s <strong> PASS ... VEC = %.3f dB</strong>\n', msg, VEC_dB);
                else
                    fprintf(2,'%s <strong> FAIL ... VEC = %.3f dB</strong>\n',  msg, VEC_dB);
                end
                
                if VEO_mV >= param.Min_VEO && VEO_mV <= param.Max_VEO
                    fprintf('%s <strong> PASS ... EH = %.3f mV</strong>\n', msg, VEO_mV);
                else
                    fprintf(2,'%s <strong> FAIL ... EH = %.3f mV</strong>\n',  msg, VEO_mV);
                end
                if COM >= param.pass_threshold
                    fprintf('%s <strong> PASS ... COM = %.3f dB</strong>\n', msg, COM);
                else
                    fprintf(2,'%s <strong> FAIL ... COM = %.3f dB</strong>\n',  msg, COM);
                end
                % begin yasuo patch 3/18/2019
                fprintf('%s DER = %.3e at COM threshold \n', msg, threshold_DER);
                % end yasuo patch
            case 'C2M'
                if VEC_dB <= param.VEC_pass_threshold
                    fprintf('%s <strong> PASS ... VEC = %.3f dB</strong>\n', msg, VEC_dB);
                else
                    fprintf(2,'%s <strong> FAIL ... VEC = %.3f dB</strong>\n',  msg, VEC_dB);
                end
                
                if VEO_mV >= param.Min_VEO && VEO_mV <= param.Max_VEO
                    fprintf('%s <strong> PASS ... EH = %.3f mV</strong>\n', msg, VEO_mV);
                else
                    fprintf(2,'%s <strong> FAIL ... EH = %.3f mV</strong>\n',  msg, VEO_mV);
                end
        end
    end
    if OP.ERL
        if ~isempty(ERL)
            if min_ERL >= param.ERL_pass_threshold
                %                             msgtext{package_testcase_i+1}=[sprintf(' PASS ... ERL = %.3f dB\n', min_ERL)];
                fprintf('%s: <strong> PASS ... ERL = %.3f dB (%.3f dB, %.3f dB)</strong>\n',msg, min_ERL, ERL );
            else
                %                             msgtext{package_testcase_i+1}=[ sprintf(' FAIL ... ERL = %.3f dB\n', min_ERL) ];
                fprintf(2,'%s:  <strong> FAIL ... ERL = %.3f dB (%.3f dB, %.3f dB)</strong>\n',msg, min_ERL, ERL) ;
            end
        end
    end
end

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

function [idx]=findbankloc(hisi,idx_st,idx_en,tap_bk,curval,bmaxg,N_bg)
% [idx]=findbankloc(hisi,idx_st,idx_en,tap_bk)
% find the location of the DFE bank
% hisi: waveform with cursor values;
% idx_st: starting index;
% idx_en: ending index ;
% tap_bk: number of taps per bank;
% bmaxg: maximum coefficient;

hisi=hisi(:);
len=idx_en-idx_st+1;
h0=abs(hisi(idx_st:idx_en));
h1=max(0,h0-bmaxg*curval);

h0n=zeros(len-tap_bk+1,1);
h1n=h0n;

for ii=1:tap_bk
    h0tmp=h0(ii:ii+len-tap_bk);
    h0n=h0n+h0tmp.^2;
    h1tmp=h1(ii:ii+len-tap_bk);
    h1n=h1n+h1tmp.^2;
end

ndiff=h0n-h1n;



idx=zeros(1,tap_bk*N_bg);
ordered_set=(1:(N_bg-1)*tap_bk+1)';
set_next_bank=0;
%Loop through each group
for k=1:N_bg
    %Sort to choose the strongest
    [~,val_sort]=sort(ndiff,'descend');
    if k==1
        %shortcut:  Choose the first 1:N_bg*tap_bk taps if they are the strongest
        if isequal(sort(val_sort(ordered_set)),ordered_set)
            idx=1:N_bg*tap_bk;
            break;
        end
    end
    if set_next_bank>0
        %when a previous bank (goodV) was found, automatically set the bank without going through the search
        new_bank=set_next_bank:set_next_bank+tap_bk-1;
        idx(tap_bk*(k-1)+1:tap_bk*k)=new_bank;
        set_next_bank=0;
        ndiff(new_bank)=0;
        bad_start=new_bank(1)-tap_bk+1;
        bad_end=new_bank(1)-1;
        if bad_end<=0
            badV=[];
        elseif bad_start>0
            badV=bad_start:bad_end;
        else
            badV=1:bad_end;
        end
        if ~isempty(badV)
            ndiff(badV)=0;
        end
        continue;
    end
    %potential bank = the strongest tap group
    new_bank=val_sort(1):val_sort(1)+tap_bk-1;
    if k==N_bg
        %Last group:  just choose the strongest
        idx(tap_bk*(k-1)+1:tap_bk*k)=new_bank;
        break;
    end
    
    do_it_again=1;
    first_time=1;
    while do_it_again
        do_it_again=0;
        %note badV:  taps smaller and less than 1 group away
        bad_start=new_bank(1)-tap_bk+1;
        bad_end=new_bank(1)-1;
        if bad_end<=0
            badV=[];
        elseif bad_start>0
            badV=bad_start:bad_end;
        else
            badV=1:bad_end;
        end
        for j=length(badV):-1:1
            if any(badV(j)-idx==0)
                badV(j)=[];
            end
        end
        %note goodV:  the tap exactly 1 tap_bk smaller
        goodV=new_bank(1)-tap_bk;
        if ~isempty(badV)
            if ~first_time
                [~,val_sort]=sort(ndiff,'descend');
            end
            first_time=0;
            checkV=[badV new_bank];

            badV_pos=zeros(1,length(badV));
            for j=1:length(badV)
                badV_pos(j)=find(badV(j)==val_sort);
            end
            
            %loop through the sorted list to find the first tap outside the group and not a member of badV
            found_goodV=0;
            for ii=1:length(val_sort)
                if val_sort(ii)==goodV
                    found_goodV=1;
                    break;
                end
                if all(val_sort(ii)-checkV~=0)
                    break;
                end
            end
            
            if ~found_goodV && min(badV_pos)<ii
                %if goodV wasn't found and bad taps occur before non group members are found
                %throw out the strongest tap and take the next strongest
                do_it_again=1;
                ndiff(new_bank(1))=0;
                %speed up:  new max_val is always val_sort(2)
                new_bank=val_sort(2):val_sort(2)+tap_bk-1;
            end
            if found_goodV
                %if goodV was found, set the next bank to goodV
                set_next_bank=goodV;
            end
        end
        
    end
    %at the end, the floating taps are set to idx
    %and ndiff has illegal values set to zero
    ndiff(new_bank)=0;
    idx(tap_bk*(k-1)+1:tap_bk*k)=new_bank;
    if ~isempty(badV)
        ndiff(badV)=0;
    end
end


idx=idx+idx_st-1;
function [tap_loc,tap_coef,hisi,b]=floatingDFE( hisi, N_b, N_bf, N_bg, N_bmax, bmaxg, curval, dfe_delta)

% hisi = postcursor isi
% N_b = number of fixed dfe taps (before floating taps begin)
% N_bf = number of floating taps per group
% N_bg = nubmber of groups
% N_bmax = max tap number that can be used for floating tap
% bmaxg = max tap strength for floating taps
% curval = value of the cursor


if nargin<8, dfe_delta=0;end


tap_coef=zeros(1,length(hisi));
b=zeros(1,length(hisi));


[tap_loc]=findbankloc(hisi,N_b+1,N_bmax,N_bf,curval,bmaxg,N_bg);

%Apply DFE to all taps
flt_curval=hisi(tap_loc);
if dfe_delta~=0
    flt_curval_q=floor(abs(flt_curval/curval)./dfe_delta) .* ...
        dfe_delta.*sign(flt_curval)*curval;
else
    flt_curval_q=hisi(tap_loc);
end
applied_coef=min(abs(flt_curval_q/curval),bmaxg).*sign(flt_curval_q);
hisi(tap_loc)= hisi(tap_loc) - curval*applied_coef;
tap_coef(tap_loc)=applied_coef;



tap_loc=sort(tap_loc,'ascend');
b(tap_loc)=bmaxg;
function [ bmax floating_tap_locations] = floating_taps_1sttest( hisi,N_b,N_bf,N_bg,N_bmax, bmaxg, COOP )
% Richard Mellitz:  04/23/2019
% hisi is the isi 1 ui/sample
% N_b number of fixed dfe taps
% N_bf number of floating taps per group
% N_bg number of floating tap groups. 1 2 or 3 right now
% N_bmax number of ui for the max reach of the floating taps
% bmaxg limit for the floating taps
% COOP = 1 co-optimize banks , -0 sequenatial optmization
%
%
% function to remove isi or add noise above bmaxg
if ~exist('COOP','var'), COOP=0;end
if iscolumn(hisi); hisi=hisi.';end
hsis_in=hisi;
% find all the reduction group taken N_bf at a time
% we are looking for the group when when remove yield the miminim isi, h, power
best_sigma=inf;best_ig1=-1;best_ig2=-1;best_ig3=-1;
% add on switch and loop for each potential group
switch N_bg
    case 0
        bmax=0;
        return
    case 1
        end1=N_bmax-N_bf;
        end2=N_b+1;
        end3=N_b+1;
    case 2
        end1=N_bmax-N_bf;
        end2=N_bmax-N_bf;
        end3=N_b+1;
    case 3
        end1=N_bmax-N_bf;
        end2=N_bmax-N_bf;
        end3=N_bmax-N_bf;
end
if COOP
    for ig1= N_b+1:end1 % now remove the 2nd group
        hcap=  hrem(hisi,ig1,N_bf,bmaxg)  ;
        % loop for 2rd group
        for ig2= N_b+1: end2
            hcap2= hrem(hcap,ig2,N_bf,bmaxg)  ;
            if N_bg < 2; hcap2 =hcap; end
            for ig3= N_b+1: end3
                hcap3= hrem(hcap2,ig3,N_bf,bmaxg)  ;
                if N_bg < 3 ; hcap3=hcap2 ; end
                sigma=norm( hcap3 );
                if sigma < best_sigma
                    best_sigma=sigma;
                    best_ig1=ig1;
                    best_ig2=ig2;
                    best_ig3=ig3;
                    best_hcap=hcap3;
                end
            end
        end
    end
else % sequentail
    for ig1= N_b+1:end1 % now remove the 1st group
        hcap=  hrem(hisi,ig1,N_bf,bmaxg)  ;
        sigma=norm( hcap );
        if sigma < best_sigma
            best_sigma=sigma;
            best_ig1=ig1;
            best_hcap=hcap;
        end
    end
    % loop for 2rd group
    hisi=best_hcap;
    for ig2= N_b+1: end2
        hcap= hrem(hisi,ig2,N_bf,bmaxg)  ;
        sigma=norm( hcap );
        if sigma < best_sigma
            best_sigma=sigma;
            best_ig2=ig2;
            best_hcap=hcap;
        end
    end
    hisi=best_hcap;
    % loop for 3rd group
    for ig3= N_b+1: end3
        hcap= hrem(hisi, ig3,N_bf,bmaxg)  ;
        sigma=norm( hcap );
        if sigma < best_sigma
            best_sigma=sigma;
            best_ig3=ig3;
            best_hcap=hcap;
        end
    end
    
end
bmax(N_b+1:N_bmax)=zeros(1,N_bmax-N_b);
switch N_bg
    case 1
        bmax(best_ig1:best_ig1+N_bf-1)=ones(1,N_bf)*bmaxg;
        floating_tap_locations= [best_ig1:best_ig1+N_bf-1];
    case 2
        bmax(best_ig1:best_ig1+N_bf-1)=ones(1,N_bf)*bmaxg;
        bmax(best_ig2:best_ig2+N_bf-1)=ones(1,N_bf)*bmaxg;
        floating_tap_locations= [best_ig1:best_ig1+N_bf-1 best_ig2:best_ig2+N_bf-1 ];
    case 3
        bmax(best_ig1:best_ig1+N_bf-1)=ones(1,N_bf)*bmaxg;
        bmax(best_ig2:best_ig2+N_bf-1)=ones(1,N_bf)*bmaxg;
        bmax(best_ig3:best_ig3+N_bf-1)=ones(1,N_bf)*bmaxg;
        floating_tap_locations= [best_ig1:best_ig1+N_bf-1 best_ig2:best_ig2+N_bf-1 best_ig3:best_ig3+N_bf-1 ];
end
floating_tap_locations=sort(floating_tap_locations);
if 0 % for code debug
    close force all
    stem(best_hcap,'disp','hcap')
    hold on
    stem(bmax,'-k','disp','bmax')
    stem(hisi,'disp','hisi')
    hold off
end

% function [ bmax floating_tap_locations] = floating_taps( hisi,   N_b,     N_bf,  N_bg,  N_bmax,       bmaxg. COO))
function [ Vfiltered, Cmod, idx ]= force( V ,param, OP , ix, C, return_V, chdata, txffe, Noise_XC)
% Vfilter is vector forced filtered sbr
% Cmod is the ffe tap co-efficient vector
% if C is passed, just process V with C else compute C
% cmx=param.rx_cmx; number of pre cursor taps
% cpx=param.rx_cps; number of post cursor taps
% V=sbr; pass pulse response
% ix the sample point in the passed pulse response
% the sample point is recomputed by optimize_fom
% idx - return floating tap location (RIM 9-19-2023)
% OP not used for now
%return_V is a flag with default value = 1.  If 0, Vfiltered is not returned
%   this allows significant speed up in optimize_fom since FFE is time consuming
%   and many combinatiFons of "Cmod" result in illegal combinations that do not need
%   Vfiltered to be calculated
% test with load('SBR_FIR_resp.mydata','-mat')
idx=[];
if nargin<4
    ix=find(V==max(V),1,'first');
end
if nargin<5
    C=[];
end
if nargin<6
    return_V=1;
end
cmx=param.RxFFE_cmx;
cpx=param.RxFFE_cpx;
% do this early on so we can reuse the old code
if param.N_bg ~=0 % must be floating taps
    cpx=param.N_bmax; % N_f in spreadsheet
end
num_taps=cmx+cpx+1;
cstep=param.RxFFE_stepz;
ndfe=param.ndfe;
spui=param.samples_per_ui;
param.current_ffegain=0;
if return_V && ~isempty(C)
    %    RIM 2-3-23 when we just want to EQ not find EQ
    Vfiltered=FFE( C , param.RxFFE_cmx,spui, V );
    Cmod=C;
    return
end
% Aling V to ix ( the sample point) and then create the sampled vector vsampled_raw
if ix < length(V)
    if isrow(V)
        if mod(ix,spui) == 0
            vsampled_raw = [V(spui+mod(ix,spui):spui:(mod(ix,spui)+spui*(floor(ix/spui)-1)))'; V(ix:spui:end)'];
        else
            vsampled_raw = [V(mod(ix,spui):spui:(mod(ix,spui)+spui*(floor(ix/spui)-1)))'; V(ix:spui:end)'];
        end

    else
        if mod(ix,spui) == 0
            vsampled_raw = [V(spui+mod(ix,spui):spui:(mod(ix,spui)+spui*(floor(ix/spui)-1))); V(ix:spui:end)];
        else
            vsampled_raw = [V(mod(ix,spui):spui:(mod(ix,spui)+spui*(floor(ix/spui)-1))); V(ix:spui:end)];
        end
    end
else
    if isrow(V)
        if mod(ix,spui) == 0%Yasou Hidaka 11/16/2018
            vsampled_raw = V(spui+mod(ix,spui):spui:end)';%Yasou Hidaka 11/16/2018
        else
            vsampled_raw = V(mod(ix,spui):spui:end)';%Yasou Hidaka 11/16/2018
        end
    else
        if mod(ix,spui) == 0%Yasou Hidaka 11/16/2018
            vsampled_raw = V(spui+mod(ix,spui):spui:end) ;
        else
            vsampled_raw = V(mod(ix,spui):spui:end) ;%Yasou Hidaka 11/16/2018
        end
    end
end
% zero pad vsampled to account for PR with short delay. RIM 10-02-2023
vsampled=[zeros(1,num_taps) vsampled_raw' zeros(1,cpx)];% pad for pre and post cursor prior to shifting

%% find the index, ivs, for the sample point but in the UI resample vector, vsampled
% ivs=find(vsampled==V(ix),1,'first');% ivs is the sample point for V
% Upen Kareti suggested fix for indexing 11/04/18
if ix < length(V)
    ivs=find(vsampled==V(ix),1,'first');% ivs is the sample point for V
else
    ivs=find(vsampled == max(vsampled),1,'first');
end


%% create VV matrix of shifted UI spaced sample of the pulse response
% only consider the VV matrix that correstonds to the FFE taps
VV=zeros(num_taps,num_taps);
for i=1:num_taps
    start_idx=ivs+i-1;
    end_idx=start_idx-num_taps+1;
    VV(:,i)=vsampled(start_idx:-1:end_idx);
end
% may want to test VV*VV' for rcond here not sure what value to use but 1e-5 is always bad
%% Apply RXFFE
if isempty(C)
    switch upper(OP.FFE_OPT_METHOD)
        case 'WIENER-HOPF'
            C= WIENER_HOPF_MMSE(vsampled ,param, OP , chdata, txffe, Noise_XC,ivs) ;
            Cmod=C(1:num_taps);
        otherwise
            % cmx+1 is the cursor or sample point
            %VV=VV(:,ivs-cmx:ivs+cpx); % only consider the VV matrix that correstonds to the FFE taps
            FV=zeros(1,num_taps); % zero the forceing vector, FV first
            FV(cmx+1)=vsampled(ivs)*10^(param.current_ffegain/20); % force the voltage at sample point
            if param.ndfe~=0 && (cpx > 0) % Yasuo Hidaka suggest fix for no postC 11/11/18
                %          FV(cmx+2)=param.bmax(1)*FV(cmx+1);   % force the post cursor to bmax if dfe exists
                FV(cmx+2)=min(param.bmax(1)*FV(cmx+1),abs(vsampled(ivs+1)))*sign(vsampled(ivs+1));
            end
            %C=((VV'*VV)^-1*VV')'*FV'; % sikve for FFE taps, C
            if diff(size(VV))==0
                %For square matrix, can solve C using simple inv(VV')*FV'
                C=VV'\FV';
            else
                %otherwise use the general solution with psuedo inverse
                %note:  this is the same as doing pinv(VV') but pinv is far slower
                % C=(inv(VV'*VV)*VV')'*FV'; % sikve for FFE taps, C
                C=(inv(VV*VV')*VV)*FV'; % sikve for FFE taps, C
            end

            Cmod=C(1:num_taps);
    end


    % added for 4.2 find floating taps with either ISI or taps
    switch lower(OP.RXFFE_FLOAT_CTL)
        case 'taps'
            [idx]=findbankloc(Cmod,param.N_tail_start,param.N_bmax,param.N_bf,Cmod(cmx+1),param.bmaxg,param.N_bg);
        otherwise
            [idx]=findbankloc(VV(cmx+1,:),param.N_tail_start,param.N_bmax,param.N_bf,Cmod(cmx+1),param.bmaxg,param.N_bg);
    end
    switch lower(OP.RXFFE_TAP_CONSTRAINT)
        case 'unity cursor'
            Cmod=Cmod/Cmod(cmx+1);
        otherwise
            Cmod=C;
    end
    if cstep ~= 0
        Cmod=floor(abs(Cmod/cstep)).*sign(Cmod)*cstep;% r250 quantize with floor ad sign(C)
    end

    if ~isempty(idx)
        idx=sort(idx);
        C1=Cmod;
        % C1(param.N_tail_start:end)=0;
        C1(param.RxFFE_cmx+param.RxFFE_cpx+2:end)=0; % zero out flosting taps
        C1(cmx+1+idx)=Cmod(cmx+1+idx);
        Cmod=C1;
    else
        % Cmod=C;
    end

    % now when ussing RxFFE floating taps need to tag stems correctly and
    % make sure DFEfloating tap code does not get exectuted

    %
else
    Cmod=C;%just us the FFE taps, C, passed for filtering
end
%%
%% filter the pulse response with the solved FFE
%  (now option to avoid this and just return Cmod for speed up)
if return_V
    Vfiltered=FFE( Cmod , param.RxFFE_cmx,spui, V );
else
    Vfiltered=[];
end
function [ILN, efit]= get_ILN(sdd21,faxis_f2)
% used for FD IL fitting
% sdd21 us a complex insertion loss
db = @(x) 20*log10(abs(x));
sdd21=squeeze(sdd21);
if  iscolumn(sdd21)
    sdd21=sdd21.';
end
fmbg=[ones(length(faxis_f2),1).*transpose(abs(sdd21))  transpose(sqrt(faxis_f2)).*transpose(abs(sdd21))  transpose(faxis_f2).*transpose(abs(sdd21)) transpose(faxis_f2.^2).*transpose(abs(sdd21)) ];
warning('off','MATLAB:nearlySingularMatrix');
LGw=transpose(abs(sdd21).*db(sdd21));
alpha = ((fmbg'*fmbg)^-1)*fmbg'*LGw;
efit=(alpha(1)+alpha(2).*sqrt(faxis_f2)+alpha(3).*faxis_f2 +faxis_f2.^2.*alpha(4)   );
ILN = db(sdd21)-efit;


function [ILN, efit, TD_ILN]= get_ILN_cmp_td(sdd21,faxis_f2,OP,param,A_T)
% Complex IL fitting
% sdd21 us a complex insertion loss
% efit and ILN are in db
% faxix_f2 needs to be at least to fb
% return reflections TD_ILN.FOM based on time domain PR fit from pulse peak
% still need to settle on voltage scaling.
% maybe db(peak/Rss

OP.interp_sparam_mag= 'trend_to_DC';
OP.interp_sparam_phase= 'interp_to_DC';
% OP.interp_sparam_mag= 'linear_trend_to_DC';
% OP.interp_sparam_phase= 'extrap_cubic_to_dc_linear_to_inf';

print_for_codereview=0;
if ~exist('A_T','var')
    A_T=1;
end

db = @(x) 20*log10(abs(x));
sdd21=squeeze(sdd21);
if  iscolumn(sdd21)
    sdd21=sdd21.';
end
fmbg=[ones(length(faxis_f2),1).*transpose(sdd21)  transpose(sqrt(faxis_f2)).*transpose(sdd21)  transpose(faxis_f2).*transpose(sdd21) transpose(faxis_f2.^2).*transpose(sdd21) ];
warning('off','MATLAB:nearlySingularMatrix');
unwraplog=log(abs(sdd21))+1i*unwrap(angle(sdd21));
LGw=transpose(sdd21.*unwraplog);
alpha = ((fmbg'*fmbg)^-1)*fmbg'*LGw;
efit_C=(alpha(1)+alpha(2).*sqrt(faxis_f2)+alpha(3).*faxis_f2 +faxis_f2.^2.*alpha(4)   );
FIT=transpose(exp(transpose(efit_C)));
efit=db(abs(FIT));
ILN = db(sdd21)-efit;
% time domain
fprintf('computing TD_ILN (dB) ...')
if exist('OP','var')
%     OP.fraction_of_F_range_start_extrap_from=.95;
    OP.impulse_response_truncation_threshold =1e-7;
    
    H_bt=Bessel_Thomson_Filter(param,faxis_f2,1);
    H_bw=Butterworth_Filter(param,faxis_f2,1);
    H_t = exp(-(pi*faxis_f2/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
    H_tw=Tukey_Window(faxis_f2,param);
    H_tw=ones(1,length(faxis_f2) );
    
   [TD_ILN.REF.FIR, ...
        TD_ILN.REF.t, ...
        TD_ILN.REF.causality_correction_dB, ...
        TD_ILN.REF.truncation_dB] = s21_to_impulse_DC(sdd21.*H_bt.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    TD_ILN.REF.PR=filter(ones(1, param.samples_per_ui), 1, TD_ILN.REF.FIR);

    [TD_ILN.FIT.FIR, ...
        TD_ILN.FIT.t, ...
        TD_ILN.FIT.causality_correction_dB, ...
        TD_ILN.FIT.truncation_dB] = s21_to_impulse_DC(FIT.*H_bt.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    TD_ILN.FIT.PR=filter(ones(1, param.samples_per_ui), 1, TD_ILN.FIT.FIR);
    ipeak=find(TD_ILN.REF.PR==max(TD_ILN.REF.PR),1,'first');
    %     NrangeUI=1000;
    %     range_end=min(min(ipeak+param.samples_per_ui*NrangeUI,length(TD_ILN.FIT.PR)-param.samples_per_ui ),length(TD_ILN.REF.PR)-param.samples_per_ui);
    range_end= min(length(TD_ILN.REF.PR), length(TD_ILN.FIT.PR));
    range=ipeak:range_end;
    TD_ILN.ILN=TD_ILN.FIT.PR(range)-TD_ILN.REF.PR(range);
    TD_ILN.t=TD_ILN.FIT.t(range);
    TD_ILN.FOM=-inf;
    TD_ILN.FOM_PDF=-inf;
    rms_fom=-inf;
    for im=1:param.samples_per_ui
        TD_ILN.FOM=max(TD_ILN.FOM, norm( TD_ILN.ILN(im:param.samples_per_ui:end)));
        [ pdf ] = get_pdf_from_sampled_signal(  TD_ILN.ILN(im:param.samples_per_ui:end), param.levels, OP.BinSize ,0);
        rms=sqrt(pdf.y*pdf.x(:).^2)*sqrt(2);
        cdf=pdf; cdf.y=cumsum(pdf.y);
        %         cursors = d_cpdf(OP.BinSize,param.a_thru*[-1 1], [1 1]/2);
        %         signal_and_isi_pdf = conv_fct(cursors, pdf);
        %         cdf=signal_and_isi_pdf; cdf.y=cumsum(signal_and_isi_pdf.y);
        if print_for_codereview  % remove once all checked out
            h=figure(190);set(gcf,'Tag','COM');
            semilogy(-cdf.x,cdf.y);
            %             xlim ([0,-cdf.x(1)])
            ylim([param.specBER 1]);title ('CDF of ILN')
            hold on
        end
        if rms>rms_fom
            rms_fom=rms;
            TD_ILN.FOM_PDF= -cdf.x(find(cdf.y >= param.specBER, 1, 'first'));
            TD_ILN.PDF=pdf;
        end
    end
    pdf_from_norm=normal_dist(TD_ILN.FOM, 7 , OP.BinSize);
    TD_ILN.SNR_ISI_FOM=db(TD_ILN.FIT.PR(ipeak)/TD_ILN.FOM);
    TD_ILN.SNR_ISI_FOM_PDF=db(TD_ILN.FIT.PR(ipeak)/TD_ILN.FOM_PDF);
    %     fprintf('%g dB\n',TD_ILN.SNR_ISI_FOM)
    fprintf('%g dB\n',TD_ILN.SNR_ISI_FOM_PDF)
    if print_for_codereview % remove once all checked out
        figure(9000);set(gcf,'Tag','COM');
        plot(TD_ILN.t,TD_ILN.ILN,'disp','td iln')
        hold on
        plot(TD_ILN.FIT.t,TD_ILN.FIT.PR,'disp','fit')
        plot(TD_ILN.REF.t,TD_ILN.REF.PR,'disp','ref')
        hold off
        fprintf('SNR ISI FOM rms = %g dB;   SNR ISI FOM PDF = %g dB\n',TD_ILN.SNR_ISI_FOM,TD_ILN.SNR_ISI_FOM_PDF)
        figure(9002);set(gcf,'Tag','COM');
        semilogy(TD_ILN.PDF.x,TD_ILN.PDF.y,'disp','actual PDF')
        hold on
        semilogy(pdf_from_norm.x,pdf_from_norm.y,'disp','PDF using Gaussian assumed PDF');
        ylim([param.specBER max([TD_ILN.PDF.y pdf_from_norm.y])]);title ('Compare actual PDF to Gaussian')
        grid on
        legend('show')
    end
end
% display('got to end of get_ILN_cmp_td')
function result = get_PSDs(result,h,cursor_i, txffe,G_DC,G_DC2,param,chdata,OP)
if 1 % force indent for doc
    num_ui=param.num_ui_RXFF_noise;
    M=param.samples_per_ui;
    L=param.levels;
    sigma_X2=(L^2-1)/(3*(L-1)^2);
    f_b=param.fb;
    T_b=1/f_b;
    delta_f = f_b/num_ui; % Units are Hz.
    fvec = (0:num_ui*M/2)*delta_f; % Single-sided frequency axis.
    result.fvec=fvec;
    SNR_TX=param.SNR_TX;
    eta_0=param.eta_0; %V^2/GHz
end
if OP.WO_TXFFE % to speed up loop find sn onlu first time ctle is case
    %% compute S_rn healey_3dj_01_2401 slide 5
    S_RN_of_f=S_RN(fvec,G_DC,G_DC2,param);
    rxn_psd=[real(S_RN_of_f(1)), S_RN_of_f(2:end-1), real(S_RN_of_f(end)), conj(S_RN_of_f(end-1:-1:2))]; % Convert single-sided frequency response to conjugate-symmetric
    rxn_psd=rxn_psd/1e9;% Units are V^2/Hz.
    rxn_rms = sqrt(sum(rxn_psd)* delta_f);
    S_rn = sum(reshape(rxn_psd, num_ui, M).');
    S_rn=S_rn(1:num_ui/2+1);
    S_rn= [real( S_rn(1)),  S_rn(2:end-1), real( S_rn(end)), conj( S_rn(end-1:-1:2))];
    %rxn_acf_samp = ifft(S_rn)*f_b; % Autocorrelation function. Note that the first term is rxn_rms^2.
    if 0 % for debug
        figure
        set(gcf, 'tag', 'COM');movegui(gcf,'northeast');
        plot(fvec(1:num_ui)/f_b,10*log10((S_rn)*1000/100) ...
            ,'disp','Srn')
        xlim([0 0.5])
        ylim([-200 -140])
        set(gcf,'defaulttextinterpreter','none')
        xlabel('Normalized Frequency')
        ylabel('PSD dBm/Hz')
        hold on
        grid on
        title('PSD')
    end
    result.S_rn=S_rn;
    result.S_rn_rms=rxn_rms;

else % find noise for item that set have tx ffe for each loop
    %% from healey_3dj_01_2401 slide 6
    % Crosstalk power spectral density
    result.S_xn=0;
    if length(chdata)~=1;
        for xchan=2:length(chdata)
            pulse_ctle=filter(ones(1,M),1,chdata(xchan).ctle_imp_response(:).');
            pulse_ctle=[ pulse_ctle(1:floor(length(pulse_ctle)/M)*M) ];
            % enable less UI for computation speed improvement
            %%
            if num_ui*M > length(pulse_ctle)
                pulse_ctle= [ pulse_ctle zeros(1,num_ui*M-length(pulse_ctle)) ];
            else
                pulse_ctle=pulse_ctle(1:num_ui*M);
            end
            cmx=find(txffe==max(txffe))-1;
            for i1=1:M
                if ~strcmp(chdata(xchan).type,'NEXT')
                    hk(xchan).k=FFE( txffe , cmx,M, pulse_ctle )'; % need to speed up here
                else
                    hk(xchan).k=pulse_ctle;
                end
            end
            for i1=1:M
                hxn(i1)=norm(hk(xchan).k(i1:M:length(hk(xchan).k)) );
            end
            iphase(xchan)=find(hxn==max(hxn));
            hk(xchan).hrn= hk(xchan).k(iphase(xchan):M:length(hk(xchan).k));
            result.hk(xchan).hrn= hk(xchan).hrn;
            hk(xchan).S_xn=sigma_X2*(abs(fft(hk(xchan).hrn))).^2/param.fb;
            result.S_xn=hk(xchan).S_xn+result.S_xn;
        end
        result.xn_rms = sqrt(sum(result.S_xn)* delta_f);
    else
        result.xn_rms=0;
        iphase=0;
    end
    %% from healey_3dj_01_2401 slide 7
    % Transmitter noise power spectral density
    htn=filter(ones(1,M),1,chdata(1).ctle_imp_response); % ctle_imp_response does not have TxFFE included
    htn=htn(mod(cursor_i,M)+1:end-mod(cursor_i,M)); % align to sample point
    htn=reshape(htn,1,[]); % make row vectors
    htn=[ htn(1:floor(length(htn)/M)*M) ];
    htn= [htn zeros(1,num_ui*M-length(htn)) ]; 
    htn=htn(1:M:end);% resample
    if num_ui>length(htn)
        hext=[htn zeros(1,num_ui-length(htn))];
    else
        hext=htn(1:num_ui);
    end
    result.S_tn=sigma_X2*10^(-SNR_TX/10)*(abs(fft(hext))).^2/param.fb; % this corresponds to +/- pi
    result.tn_rms = sqrt(sum(result.S_tn)* delta_f); 
    %% from healey_3dj_01_2401 slide 8
    % Power spectral density of noise due to jitter
    %% Eq. 93A-28 %%
    sampling_offset = mod(cursor_i, M);
    %ensure we can take early sample
    if sampling_offset<=1
        sampling_offset=sampling_offset+M;
    end
    if (OP.LIMIT_JITTER_CONTRIB_TO_DFE_SPAN)
        cursors_early_sample = h(cursor_i-1+M*(-1:param.ndfe));
        cursors_late_sample = h(cursor_i+1+M*(-1:param.ndfe));
    else
        cursors_early_sample = h(sampling_offset-1:M:end);
        cursors_late_sample = h(sampling_offset+1:M:end);
    end
    % ensure lengths are equal
    cursors_early_sample = cursors_early_sample(1:length(cursors_late_sample));
    h_J = (cursors_late_sample-cursors_early_sample)/2*M;
    h_J=reshape(h_J,1,[]); % make row vectors
    if num_ui>length(h_J)
        h_J=[h_J zeros(1,num_ui-length(h_J))];
    else
        h_J=h_J(1:num_ui);
    end
    result.iphase=iphase;
    result.S_jn=sigma_X2*(param.A_DD^2+param.sigma_RJ^2)*(abs(fft(h_J))).^2/param.fb; % this corresponds to +/- pi
    result.jn_rms = sqrt(sum(result.S_jn)* delta_f);
    result.S_n=result.S_rn+ result.S_tn+ result.S_xn+ result.S_jn;
end
function result=get_PulseR(ir,param,cb_step,ZT)
%ir = impulse response
%t_base=time array with equal time steps
%samp_UI = number of samples per UI for ir

% t for debug
t=(1/param.fb)/param.samples_per_ui*(0:length(ir)-1);

if cb_step
    Ag=1;
    dt=1/param.fb/param.samples_per_ui;
    edge_time=param.TR_TDR*1e-9;
    fedge=1/edge_time;
    tedge=0:dt:edge_time*2;
    %
    edge=Ag*(2*cos(2*pi*(tedge)*fedge/16-pi/4).^2-1);
    drive_pulse=[edge ones(1,param.samples_per_ui)];
    %pulse=filter(UI_ones,1,ir);
    % t for debug
    t=(1/param.fb)/param.samples_per_ui*(0:length(ir)-1);
    
    pulse=filter(drive_pulse,1,ir);
else
    pulse=filter( ones(1,param.samples_per_ui),1,ir);
end
PDR_response=(1+pulse)./(1-pulse).*ZT*2;
result.PDR=PDR_response;
result.pulse=pulse;



function [ FIR t] =get_RAW_FIR(H,f,OP,param)
H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*f./(0.75*param.fb));
if ~iscolumn(H), H=H.';end
if ~iscolumn(H_r), H_r=H_r.';end
H=H(:).*H_r;
[FIR, t, ~,~] = s21_to_impulse_DC(H ,f, param.sample_dt, OP) ;
% SBR=filter(ones(1, param.samples_per_ui), 1, FIR);

function RILN_TD_struct= get_RILN_cmp_td(sdd21,RIL_struct,faxis_f2,OP,param,A_T)
% Complex reflection and re-reflection noise using the concept of zero'ing
% out of reflections
% sdd21 us a complex insertion loss
% RIL_struct is the output of capture_RIL_RILN()
% faxix_f2 needs to be at least to fb
% return reflections RILN_TD_struct.FOM based on time domain PR fit from pulse peak
% still need to settle on voltage scaling.
% maybe db(peak/Rss
db = @(x) 20*log10(abs(x));
fprintf('computing TD_RILN (dB) ...');

OP.interp_sparam_mag= 'trend_to_DC';
OP.interp_sparam_phase= 'interp_to_DC';
% OP.interp_sparam_mag= 'linear_trend_to_DC';
% OP.interp_sparam_phase= 'extrap_cubic_to_dc_linear_to_inf';

sdd21=squeeze(sdd21);
if iscolumn(sdd21)
    sdd21=sdd21.';
end
RIL=squeeze(RIL_struct.RIL);
if  iscolumn(RIL)
    RIL=RIL.';
end
rho_port1=squeeze(RIL_struct.rho_port1);
if  iscolumn(rho_port1)
    rho_port1=rho_port1.';
end
rho_port2=squeeze(RIL_struct.rho_port2);
if  iscolumn(rho_port2)
    rho_port2=rho_port2.';
end
RIL_f=squeeze(RIL_struct.freq);
if  iscolumn(RIL_f)
    RIL_f=RIL_f.';
end

%---start. Calculate the reflection and re-reflection noise
number_of_echos= 1e3;
fmin= 1e9;%<-------------
port2_reflection_rereflection_noise= zeros(1, length(RIL));
port1_reflection_rereflection_noise= zeros(1, length(RIL));
for m= 1:number_of_echos
    port2_reflection_rereflection_noise= port2_reflection_rereflection_noise+ abs(RIL).*(RIL.^(2*m)).*(rho_port1.^m).*(rho_port2.^m).*(1+rho_port1).*(1+rho_port2);
    port1_reflection_rereflection_noise= port1_reflection_rereflection_noise+ abs(RIL).*(RIL.^(2*m-1)).*(rho_port1.^(m-1)).*(rho_port2.^m).*(1+rho_port1).*(1+rho_port1);
end

%-----start. In the case of reflections, observed is bad TD conversion and hence removing data before 1GHz
fmin_idx= find(RIL_f>= fmin, 1, 'first');
port2_reflection_rereflection_noise= port2_reflection_rereflection_noise(fmin_idx:end);
port1_reflection_rereflection_noise= port1_reflection_rereflection_noise(fmin_idx:end);
f_reflection_rereflection_noise= RIL_f(fmin_idx:end);
%-----end. In the case of reflections, observed is bad TD conversion and hence removing data before 1GHz

% clear RIL RIL_f rho_port1 rho_port2
% clear fmin m
%---end. Calculate the reflection and re-reflection noise

fmbg=[ones(length(faxis_f2),1).*transpose(sdd21)  transpose(sqrt(faxis_f2)).*transpose(sdd21)  transpose(faxis_f2).*transpose(sdd21) transpose(faxis_f2.^2).*transpose(sdd21) ];
warning('off','MATLAB:nearlySingularMatrix');
unwraplog=log(abs(sdd21))+1i*unwrap(angle(sdd21));
LGw=transpose(sdd21.*unwraplog);
alpha = ((fmbg'*fmbg)^-1)*fmbg'*LGw;
efit_C=(alpha(1)+alpha(2).*sqrt(faxis_f2)+alpha(3).*faxis_f2 +faxis_f2.^2.*alpha(4)   );
FIT=transpose(exp(transpose(efit_C)));
efit=db(abs(FIT));
ILN = db(sdd21)-efit;


OP.impulse_response_truncation_threshold =1e-7;

print_for_codereview=0;
if exist('OP','var')
    H_bt=Bessel_Thomson_Filter(param,faxis_f2,1);
    H_bw=Butterworth_Filter(param,faxis_f2,1);
    H_t = exp(-(pi*faxis_f2/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
    H_tw=Tukey_Window(faxis_f2,param);
    H_tw=ones(1,length(faxis_f2) );
    [RILN_TD_struct.REF.FIR, ...
        RILN_TD_struct.REF.t, ...
        RILN_TD_struct.REF.causality_correction_dB, ...
        RILN_TD_struct.REF.truncation_dB] = s21_to_impulse_DC(sdd21.*H_bw.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    RILN_TD_struct.REF.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.REF.FIR);
    
        
    [RILN_TD_struct.FIT.FIR, ...
        RILN_TD_struct.FIT.t, ...
        RILN_TD_struct.FIT.causality_correction_dB, ...
        RILN_TD_struct.FIT.truncation_dB] = s21_to_impulse_DC(FIT.*H_bw.*H_t.*H_tw ,faxis_f2, param.sample_dt, OP) ;
    RILN_TD_struct.FIT.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.FIT.FIR);
    
    
    H_bt=Bessel_Thomson_Filter(param,RIL_f,1);
    H_bw=Butterworth_Filter(param,RIL_f,1);
    H_t = exp(-(pi*RIL_f/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
    H_tw=Tukey_Window(RIL_f,param);
    H_tw=ones(1,length(RIL_f) );
    [RILN_TD_struct.RIL.FIR, ...
        RILN_TD_struct.RIL.t, ...
        RILN_TD_struct.RIL.causality_correction_dB, ...
        RILN_TD_struct.RIL.truncation_dB] = s21_to_impulse_DC(RIL.*H_bw.*H_t.*H_tw ,RIL_f, param.sample_dt, OP) ;
    RILN_TD_struct.RIL.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.RIL.FIR);
    
    
    %---start. Calculate the channel delay
    try
    [delay_sec, delay_idx]= calculate_delay_CausalityEnforcement(faxis_f2, sdd21, param, OP);
    catch
    end
    port2_reflection_rereflection_noise= port2_reflection_rereflection_noise.*exp(-1j*2*pi*f_reflection_rereflection_noise*delay_sec);
    clear delay_sec delay_idx
    %---end. Calculate the channel delay


    
    H_bt=Bessel_Thomson_Filter(param,f_reflection_rereflection_noise,1);
    H_bw=Butterworth_Filter(param,f_reflection_rereflection_noise,1);
    H_t = exp(-(pi*f_reflection_rereflection_noise/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
    H_tw=Tukey_Window(f_reflection_rereflection_noise,param);
    H_tw=ones(1,length(f_reflection_rereflection_noise) );
    [RILN_TD_struct.REF_noise.FIR, ...
        RILN_TD_struct.REF_noise.t, ...
        RILN_TD_struct.REF_noise.causality_correction_dB, ...
        RILN_TD_struct.REF_noise.truncation_dB] = s21_to_impulse_DC(port2_reflection_rereflection_noise.*H_bw.*H_t.*H_tw ,f_reflection_rereflection_noise, param.sample_dt, OP) ;
    RILN_TD_struct.REF_noise.PR=filter(ones(1, param.samples_per_ui), 1, RILN_TD_struct.REF_noise.FIR);
    
    ipeak=find(RILN_TD_struct.REF.PR==max(RILN_TD_struct.REF.PR),1,'first');
    %     NrangeUI=1000;
    %     range_end=min(min(ipeak+param.samples_per_ui*NrangeUI,length(RILN_TD_struct.FIT.PR)-param.samples_per_ui ),length(RILN_TD_struct.REF.PR)-param.samples_per_ui);
    range_end= min(length(RILN_TD_struct.REF.PR), length(RILN_TD_struct.REF_noise.PR));
    range=ipeak:range_end;
    RILN_TD_struct.ILN=RILN_TD_struct.REF_noise.PR(range);
    RILN_TD_struct.t=RILN_TD_struct.REF_noise.t(range);
    RILN_TD_struct.FOM=-inf;
    RILN_TD_struct.FOM_PDF=-inf;
    rms_fom=-inf;
    for im=1:param.samples_per_ui
        RILN_TD_struct.FOM=max(RILN_TD_struct.FOM, norm( RILN_TD_struct.ILN(im:param.samples_per_ui:end)));
        [ pdf ] = get_pdf_from_sampled_signal(  RILN_TD_struct.ILN(im:param.samples_per_ui:end), param.levels, OP.BinSize ,0);
        rms=sqrt(pdf.y*pdf.x(:).^2)*sqrt(2);
        cdf=pdf; cdf.y=cumsum(pdf.y);
        %         cursors = d_cpdf(OP.BinSize,param.a_thru*[-1 1], [1 1]/2);
        %         signal_and_isi_pdf = conv_fct(cursors, pdf);
        %         cdf=signal_and_isi_pdf; cdf.y=cumsum(signal_and_isi_pdf.y);
        if print_for_codereview  % remove once all checked out
            h=figure(190);set(gcf,'Tag','COM');
            semilogy(-cdf.x,cdf.y);
            %             xlim ([0,-cdf.x(1)])
            ylim([param.specBER 1]);title ('CDF of ILN')
            hold on
        end
        if rms>rms_fom
            rms_fom=rms;
            RILN_TD_struct.FOM_PDF= -cdf.x(find(cdf.y >= param.specBER, 1, 'first'));
            RILN_TD_struct.PDF=pdf;
        end
    end
    pdf_from_norm=normal_dist(RILN_TD_struct.FOM, 7 , OP.BinSize);
    RILN_TD_struct.SNR_ISI_FOM=db(RILN_TD_struct.FIT.PR(ipeak)/RILN_TD_struct.FOM);
    RILN_TD_struct.SNR_ISI_FOM_PDF=db(RILN_TD_struct.FIT.PR(ipeak)/RILN_TD_struct.FOM_PDF);
    %     fprintf('%g dB\n',RILN_TD_struct.SNR_ISI_FOM)
    fprintf('%g dB\n',RILN_TD_struct.SNR_ISI_FOM_PDF)
    if print_for_codereview % remove once all checked out
        figure(9000);set(gcf,'Tag','COM');
        plot(RILN_TD_struct.REF.t,RILN_TD_struct.REF.PR,'disp','ref')
        hold on
        plot(RILN_TD_struct.FIT.t,RILN_TD_struct.FIT.PR,'disp','fit')
        plot(RILN_TD_struct.RIL.t,RILN_TD_struct.RIL.PR,'disp','RILN')
        yyaxis right
        plot(RILN_TD_struct.t,RILN_TD_struct.ILN,'disp','td RILN')
        hold off
        fprintf('SNR ISI FOM rms = %g dB;   SNR ISI FOM PDF = %g dB\n',RILN_TD_struct.SNR_ISI_FOM,RILN_TD_struct.SNR_ISI_FOM_PDF)
        figure(9002);set(gcf,'Tag','COM');
        semilogy(RILN_TD_struct.PDF.x,RILN_TD_struct.PDF.y,'disp','actual PDF')
        hold on
        semilogy(pdf_from_norm.x,pdf_from_norm.y,'disp','PDF using Gaussian assumed PDF');
        ylim([param.specBER max([RILN_TD_struct.PDF.y pdf_from_norm.y])]);title ('Compare actual PDF to Gaussian')
        grid on
        legend('show')
    end
end
function result=get_StepR(ir,param,cb_step,ZT)
%ir = impulse response
%t_base=time array with equal time steps
%samp_UI = number of samples per UI for ir
%   result.SBR
% t for debug
t=(1/param.fb)/param.samples_per_ui*(0:length(ir)-1);

if cb_step
    Ag=1;
    dt=1/param.fb/param.samples_per_ui;
    edge_time=param.TR_TDR*1e-9;
    fedge=1/edge_time;
    tedge=0:dt:edge_time*2;
    %
    edge=Ag*(2*cos(2*pi*(tedge)*fedge/16-pi/4).^2-1);
    drive_pulse=[edge ones(1,param.samples_per_ui)];
    %pulse=filter(UI_ones,1,ir);
    
    pulse=filter(drive_pulse,1,ir);
else
    pulse=cumsum(ir);
end
TDR_response=(1+pulse)./(1-pulse)*ZT*2;
result.ZSR=TDR_response;
result.pulse=pulse;


function TDR_results = get_TDR(sdd, OP, param,ZT,np)
% sdd is differential s-parameters structure (2 port assumed)
% input parameter structure for s parameters sdd--> sdd.Impedance, sdd.Frequencies, sdd.Parameters, sdd.NumPorts
% TDR_results.delay             pre t=0 delay for TDR... help with time domain responce quaility
% TDR_results.tdr               the TDR responce (ohms vs  TDR_results.t
% TDR_results.t                 starting at t=0
% TDR_results.tx_filter         transmitter filter vs TDR_results.f
% TDR_results.Rx_filter         receiver filter vs TDR_results.f
% TDR_results.f                 frequency for filter and s parameters
% TDR_results.ptdr_RL           reflection waveform from the pulse
% TDR_results.WC_ptdr_samples_t worst case time sample of the reflection pulse
% TDR_results.WC_ptdr_samples   worst case reflection samples of the reflection pulse
% TDR_results.ERL               reported effective return loss
%
db = @(x) 20*log10(abs(x));
rms =@(x) norm(x)/sqrt(length(x));
if isfield(OP,'TDR_duration')
    TDR_duration=OP.TDR_duration; % approximate transit time multipler
else
    TDR_duration=5;
end
if ~isfield(OP,'DISPLAY_WINDOW')
    OP.DISPLAY_WINDOW=1; % approximate transit time multipler
end
f=sdd.Frequencies;
TDR_results.f=f;
% OP.Zt_adj=2;
if param.FLAG.S2P == 0
    
    % re-normalize reference of s-parameterss: this seems correct for a s4p input file
    TDR_RL =@(Zin,Zout,s11,s12,s21,s22)(Zin.^2.*s11+Zin.^2.*s22+Zout.^2.*s11+Zout.^2.*s22+Zin.^2-Zout.^2+Zin.*Zout.*s11.*2.0-Zin.*Zout.*s22.*2.0+Zin.^2.*s11.*s22-Zin.^2.*s12.*s21-Zout.^2.*s11.*s22+Zout.^2.*s12.*s21)./(Zin.*Zout.*2.0+Zin.^2.*s11+Zin.^2.*s22-Zout.^2.*s11-Zout.^2.*s22+Zin.^2+Zout.^2+Zin.^2.*s11.*s22-Zin.^2.*s12.*s21+Zout.^2.*s11.*s22-Zout.^2.*s12.*s21-Zin.*Zout.*s11.*s22.*2.0+Zin.*Zout.*s12.*s21.*2.0);
    
    if param.RL_sel==1, other_port=2;end
    if param.RL_sel==2, other_port=1;end
    for i = 1:length(sdd.Frequencies)
        if size(sdd.Parameters,2) ==1 % for s2p files
            RL(i)=TDR_RL(sdd.Impedance,2*ZT,sdd.Parameters(param.RL_sel,param.RL_sel,i) ,1, 1 ,sdd.Parameters(param.RL_sel,param.RL_sel,i) );
        else
            RL(i)=TDR_RL(sdd.Impedance,2*ZT,sdd.Parameters(param.RL_sel,param.RL_sel,i) ,sdd.Parameters( 1,2 ,i), sdd.Parameters( 2,1 ,i),sdd.Parameters( other_port,other_port ,i) );
        end
    end
    % elseif OP.Zt_adj ==2 % only adjust z_t drive impedance
    %     RL=squeeze((sdd.Parameters(param.RL_sel,param.RL_sel,:)));
    %     Z_t=ZT;
    %     zref=sdd.Impedance/2;
    %     if Z_t > zref
    %         radjust=  (zref-Z_t);
    %         S11adjust= radjust./(radjust + 2*zref);
    %         RL=RL  +S11adjust;
    %     elseif Z_t < zref
    %         rpad=-Z_t*zref/(Z_t-zref);
    %         S11adjust=zref/(rpad*(zref/rpad + 2));
    %         RL=RL   + S11adjust;
    %     else
    %         RL=RL;
    %     end
else
    for i = 1:length(sdd.Frequencies)
        rho= (2*ZT - sdd.Impedance) / (2*ZT + sdd.Impedance);
        interim =  sqrt(1-abs(rho)^2) * (1 - rho) / abs(1-rho);
        RL(i) = interim \ (sdd.Parameters(param.RL_sel,param.RL_sel,i) - rho) / ...
            (1 - rho *sdd.Parameters(param.RL_sel,param.RL_sel,i) ) * interim;
    end
end

% end
RL=squeeze(RL);
f9=f/1e9;
tr=param.TR_TDR;
TDR_results.delay=500e-12 ;
% determine max time from thue
% if sdd.NumPorts==1
%     try
%         maxtime=OP.N*param.ui;
%     catch
%         maxtime=2e-9;
%     end
%     pix=1;
% else
%     [ fir4del, tu] =get_RAW_FSIR(squeeze(sdd.Parameters(2,1,:)),f,OP,param);
%     pix=find(fir4del==max(fir4del),1);
%     maxtime=tu(pix)*TDR_duration+TDR_results.delay;
%     if maxtime > tu(end); maxtime=tu(end);end
% endS

try
    maxtime=OP.N*param.ui;
catch
    maxtime=2e-9;
end
if OP.N==0
    if sdd.NumPorts==1
        fprintf('<strong> Warning for s2p files N must not be zero<\strong> ');
    else
        [ fir4del, tu] =get_RAW_FIR(squeeze(sdd.Parameters(2,1,:)),f,OP,param);
        pix=find(fir4del==max(fir4del),1);
        maxtime=tu(pix)*TDR_duration+TDR_results.delay;
        if maxtime > tu(end); maxtime=tu(end);end
    end
end


% add delay 500 ps for TDR and 3 times Gaussnan transtion time
% (makes gausian edge somewhat causal)
H_t = exp( -2*(pi*f9*(tr)/1.6832).^2 ).*exp(-(1j)*2*pi*f9*TDR_results.delay/1e-9).*exp(-1j*2*pi*f9*tr*3);
if ~isfield(OP,'cb_Guassian')
    Use_gaussian=1;
else
    Use_gaussian=OP.cb_Guassian;
end
if Use_gaussian
    if iscolumn(H_t), H_t=H_t.'; end
    RLf=RL(:).'.*H_t;
else % add extra 3x tr delay for causality
    RLf=RL.*exp(-(1j)*2*pi*f9*TDR_results.delay/1e-9);
end

%Bessesl-Thomson turned off here (3rd input=0)
H_bt=Bessel_Thomson_Filter(param,f,0);

if isfield(OP,'TDR_Butterworth')
    H_bw=Butterworth_Filter(param,f,OP.TDR_Butterworth);
else
    H_bw=ones(1,length(f));
end


if param.Tukey_Window ~= 0
    H_tw= Tukey_Window(f,param);
else
    H_tw=ones(1,length(f));
end


if iscolumn(H_tw), H_tw=H_tw.';end
if iscolumn(H_bt), H_bt=H_bt.';end
if iscolumn(H_bw), H_bw=H_bw.';end
if iscolumn(RLf), RLf=RLf.';end

TDR_results.Rx_filter=H_bt.*H_bw.*H_tw;
RLf=RLf.*TDR_results.Rx_filter;
TDR_results.tx_filter=H_t;


[IR, t, causality_correction_dB, truncation_dB] = ...
    s21_to_impulse_DC(RLf,  sdd.Frequencies(:), param.sample_dt,OP);


%
% param.tfx =4.2e-10; % need to put in xls file and comment this out
tfx=param.tfx(np); % use fixture delay for port (np)

%  IR(abs(IR)<=OP.impulse_response_truncation_threshold)=0; % noise filter

t = t-TDR_results.delay;
tend=find(t>=maxtime+tfx,1); % n starts at tfx
if isempty(tend), tend=length(t); end
IR=IR(1:tend);
t=t(1:tend);
if isempty(tend), tend=length(t); end
tstart=find(t>=tr*1e-9,1); % account for Gaussian precursor
if isempty(tstart), tstart=1;end
if isempty(tend) || tstart >= tend
    if isempty(tend) || tstart >= tend
        %         warndlg('TDR compuation not valid, try decreasing truncation tolerance, increasing samples, or adding a transmisson line','WrnTDR');
    end
    tend=length(t);
    tstart=1;
end
OP.cb_step=0; % step is a basically a cos^2 form -pi/4 to 0. not activated.
ch=get_StepR(IR(tstart:tend),param,OP.cb_step,ZT);
TDR_results.tdr= ch.ZSR;
TDR_results.t = t(tstart:tend);

PTDR=get_PulseR(IR(tstart:tend),param,OP.cb_step,ZT);
if OP.TDR || OP.PTDR % determin average impededance with
    try
        tfstart=find(t>=3*tr*1e-9,1);
        x=squeeze(TDR_results.t(tfstart:end));TDR_results.x=TDR_results.tdr(:);
        y=squeeze(TDR_results.tdr(tfstart:end));TDR_results.y=TDR_results.t(:);
        w= exp(-(x-x(1))/OP.T_k ) ; % weighting function
        TDR_results.avgZport=mean(y.*w.')/mean(w.');
    catch
        TDR_results.avgZport=0;
        fit=zeros(1,1);
        p=[0 0 0 0  ];
    end
    TDR_results.RL=RL;
end
if OP.PTDR
    %     param.N_bx=param.ndfe;
    RL_equiv=-inf;
    L=param.levels;
    BinSize=OP.BinSize;
    %     param.specBER=1e-5;
    if OP.DISPLAY_WINDOW
        hwaitbar=waitbar(0);
    else
        fprintf('Worst ERL searching');
    end
    % adjust PTDR for NDFE
    % ---------------------- 2.7 code
    %     ntx=find(TDR_results.t >= tfx,1,'first');
    %     %     gatestartt=TDR_results.t(ntx);
    %     %     gatestartV=PTDR.pulse(ntx);
    %     ndfex=find(TDR_results.t > (param.N_bx+1)*param.ui+tfx,1,'first');
    %     tk=param.ui*1*(param.N_bx+1)+tfx;
    % -------------------
    % [ahealey 09/06/2019] Need to account for the 3*TR_TDR delay included in the rise
    % time filter.
    % ntx=find(TDR_results.t >= tfx,1,'first');
    ntx=find(TDR_results.t >= tfx+3*tr*1e-9,1,'first');
    %     gatestartt=TDR_results.t(ntx);
    %     gatestartV=PTDR.pulse(ntx);
    % ndfex=find(TDR_results.t > (param.N_bx+1)*param.ui+tfx,1,'first');
    % tk=param.ui*1*(param.N_bx+1)+tfx;
    ndfex=find(TDR_results.t > (param.N_bx+1)*param.ui+tfx+3*tr*1e-9,1,'first');
    tk=param.ui*1*(param.N_bx+1)+tfx+3*tr*1e-9;
    % [ahealey] End of modifications.
    if  isempty(ndfex), ndfex=length(TDR_results.t); end
    PTDR.pulse_orig=PTDR.pulse;
    
    switch param.Grr
        case 0 % pre .3cd release
            fctrx(1:length(PTDR.pulse_orig))=(1+param.rho_x)*param.rho_x;
        case 1 % .3cd release
            fctrx(1:length(PTDR.pulse_orig))=1;
        case 2 % .3ck working
            fctrx(1:length(PTDR.pulse_orig))=1;
    end
    Gloss(1:length(TDR_results.t))=1;
    Grr(1:length(TDR_results.t))=1;
    fctrx(1:ntx)=0; % moved out of loop for beta x and n_bx
    
    for ii=ntx:ndfex
        % adjust for near end loss
        if param.N_bx>0 && param.beta_x~=0;
            Gloss(ii)= 10.^(param.beta_x*(TDR_results.t(ii)-tk)/20);
        else
            Gloss(ii)=1;
        end
        % ---------------------- 2.7 code
        %         x=(TDR_results.t(ii)-tfx)/param.ui;
        % ----------------------
        % [ahealey9/06/2019] Need to account for the 3*TR_TDR delay included in the
        % rise time filter.
        % x=(TDR_results.t(ii)-tfx)/param.ui;
        x=(TDR_results.t(ii)-tfx-3*tr*1e-9)/param.ui;
        % determine how much of the return loss to use base on expected
        % missing reflections
        switch param.Grr
            case 0 % pre .3cd release
                Grr(ii)= (1+param.rho_x)*param.rho_x*exp(-(x-1*param.N_bx-1).^2/(1+param.N_bx)^2);
            case 1 % .3cd release
                Grr(ii)= (1+param.rho_x)*param.rho_x*exp(-(x-1*param.N_bx-1).^2/(1+param.N_bx)^2);
            case 2 % .3ck working
                Grr(ii)= param.rho_x ;
        end
        fctrx(ii)=Gloss(ii).*Grr(ii);
    end
    
    if isrow(fctrx), fctrx=fctrx(:);end
    PTDR.pulse=PTDR.pulse.*fctrx;
    if 0
        figure(10101+param.RL_sel);set(gcf,'Tag','COM');
        s1=subplot(2,1,1);
        plot((TDR_results.t(ntx:end)- TDR_results.t(ntx))/param.ui,fctrx(ntx:end),'disp','G_l_o_s_s*G_r_r');
        hold on
        plot((TDR_results.t(ntx:end)- TDR_results.t(ntx))/param.ui,Gloss(ntx:end),'disp','G_l_o_s_s');
        plot((TDR_results.t(ntx:end)- TDR_results.t(ntx))/param.ui,Grr(ntx:end),'disp','G_r_r');
        grid on
        ylim([ 0 1.2])
        s2=subplot(2,1,2);
        plot((TDR_results.t(ntx:end)- TDR_results.t(ntx))/param.ui,PTDR.pulse(ntx:end)./fctrx(ntx:end),'disp','PTDR');
        grid on
        linkaxes([s1,s2],'x')
        xlabel 'UI'
        xlim ([ 1 200])
    end
    
    FAST_NOISE_CONV=0;
    ERLRMS=rms(PTDR.pulse);
    for ki=1:param.samples_per_ui
        progress = ki/param.samples_per_ui;
        if OP.DISPLAY_WINDOW
            waitbar(progress, hwaitbar, 'ERL computing'); figure(hwaitbar); drawnow;
        else
            if ~mod(progress*100,1), fprintf('%i%% ', progress*100 );end
        end
        tps=PTDR.pulse(ki:param.samples_per_ui:end);
        if OP.RL_norm_test
            rl_fom=(norm(tps));
        else
            testpdf=get_pdf_from_sampled_signal( tps,L, BinSize*10, FAST_NOISE_CONV );
            cdf_test=cumsum(testpdf.y);
            rl_test=(-testpdf.x(find(cdf_test>=param.specBER,1,'first')));
            rl_fom=rl_test;
        end
        if rl_fom > RL_equiv
            RL_equiv=rl_fom;
            best_ki=ki;
        end
        if ~OP.RL_norm_test
            best_erl=rl_test;
            best_pdf=testpdf;
            best_cdf=cdf_test;
        end
        
    end
    if OP.RL_norm_test
        tps=PTDR.pulse(best_ki:param.samples_per_ui:end);
        testpdf=get_pdf_from_sampled_signal( tps,L, BinSize*10, FAST_NOISE_CONV );
        cdf_test=cumsum(testpdf.y);
        best_erl=(-testpdf.x(find(cdf_test>=param.specBER,1,'first')));
    end
    
    fprintf('\n');
    try
        close(hwaitbar)
    catch
    end
    if ~exist('best_ki','var'),best_ki=1;end
    TDR_results.ptdr_RL=PTDR.pulse;  % reflection waveform from the pulse
    TDR_results.WC_ptdr_samples_t=TDR_results.t(best_ki:param.samples_per_ui:end);
    TDR_results.WC_ptdr_samples=PTDR.pulse(best_ki:param.samples_per_ui:end);
    TDR_results.ERL=-db(best_erl);
    TDR_results.ERLRMS=-db(ERLRMS);
    
end


% end get TDR
%%




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
function [chdata, param] = get_s4p_files(param, OP, num_fext, num_next, file_list)
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
        set(h,'Tag','COM') % RIM 06-13-2022 ... tak msg box for closing later
    end
    if  OP.ERL == 2
        dir=fullfile(filepath, '*.s2p');
        [basename,filepath]=uigetfile(dir,'input RL measurement .s2p ');
        if filepath == 0
            error('No RL measurement file')
        end
    else
        dir=fullfile(filepath, '*.s4p');
        [basename,filepath]=uigetfile(dir,'input thru channel response .s4p ');
        if filepath == 0
            error('No Thru file')
        end
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
            set(h,'Tag','COM') % RIM 06-13-2022 ... tak msg box for closing later
        end
        if  param.tfx(1) == -1 && OP.ERL_ONLY && OP.ERL == 2;
            dir=fullfile(filepath, '*.s4p');
            [basename,filepath]=uigetfile(dir,'input noise channel response .s4p');
            if filepath==0
                error('Not enough NEXT files')
            end
        else
            dir=fullfile(filepath, '*.s4p');
            if OP.RX_CALIBRATION == 1 || OP.FIXTURE_CALIBRATION
                [basename,filepath]=uigetfile(dir,'input noise channel response .s4p');
            else
                [basename,filepath]=uigetfile(dir,['input fext channel response .s4p #', num2str(nxi-kxi)]);
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
        dir=fullfile(filepath, '*.s4p');
        [basename,filepath]=uigetfile(dir,['input next channel response .s4p ', num2str(nxi-kxi)]);
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

function [sigma_N] = get_sigma_eta_ACCM_noise(chdata,param,H_sy,H_r,H_ctf)
% 11-25-2020 correct integratation limits, should be infinity or highest specified frequency
% H_sy - PSD for Tx power delivery, not normally used and set to a vector 1's
% H_r - receiver filter, Butterworth
% H_ctf - total gain of CTLE and low freq filtering
% H_dc - the common mode channel gain
% param.eta_0 -input referred Rx noise
% param.AC_CM_RMS - AC CM source before Tx series source resistor.
% param.ACCM_MAX_Freq - Max frequency for ACCM source intergration
%% Equation 93A-35 - independent of FFE setting %%
sigma_N1 = sqrt(param.eta_0*sum( abs(H_sy(2:end) .* H_r(2:end) .* H_ctf(2:end) ).^2 .* diff(chdata(1).faxis)/1e9));% eta_0 is V^2/Ghz i.e. /1
if sum(param.AC_CM_RMS) ~= 0
     sigma_ACCM=0;
    f_int= chdata(1).faxis( chdata(1).faxis<=param.ACCM_MAX_Freq );
    for i=1:length(chdata)       
        H_dc=abs(squeeze(chdata(i).sdc21));
        sigma_ACCM_acc= sqrt( 2*param.AC_CM_RMS_TX^2.*sum( abs(H_sy(2:length(f_int)) .* H_r(2:length(f_int)) .* H_ctf(2:length(f_int)).*H_dc(2:length(f_int)) ).^2 .* diff(f_int))/f_int(end) );
        sigma_ACCM=norm([sigma_ACCM_acc,sigma_ACCM]);
    end
    sigma_N=norm([sigma_N1,sigma_ACCM]);
else
    sigma_N=sigma_N1;
end
%%
function [ sigma_NE , sigma_HP] = get_sigma_noise( H_ctf, param, chdata, sigma_bn )
% for Rx calibratrion only
% the FEXT channel for calibration basically a DC connection unlike normal
% FEXT channels which are nearly open at DC channels
H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(2).faxis/(param.f_r*param.fb));
idxfbby2=find( chdata(2).faxis(:) >= param.fb/2, 1);
if size(chdata,2) >= 2
    Hnoise_channel=chdata(2).sdd21;% rx package is already included tx is not
else
    Hnoise_channel=1;
end
f=chdata(2).faxis;
f_hp=param.f_hp;
if f_hp ~=0 % param.f_hp is a key indicating that Tx bbn is used as in clause 162 
    H_hp=(-1j*f./f_hp)./(1+1j*f./f_hp);
else
    H_hp=ones(1,length(f));
end
%% Equation 93A-47 or 162-12 
H_np=Hnoise_channel.*H_ctf.*H_r.*H_hp;

%% Equation 93A-48 or 162-14%%
sigma_NE = sigma_bn*sqrt(mean(abs(H_np(1:idxfbby2).^2)));
sigma_HP = sigma_bn*(mean(abs(H_hp(1:idxfbby2).^2)));
function [sigma_XT, sigma_FEXT, sigma_NEXT] = get_xtlk_noise( upsampled_txffe, type, param, chdata, phase_memory,C )
% Modified not to double count crosstalk: John Ewen 13/12/2018
% function sigma_XT = get_xtlk_noise( upsampled_txffe, type, param, chdata,ctle_indx,clow_indx, C,cursor_i)
index_f2=find(chdata(1).faxis(:)>param.fb,1,'first');
if isempty(index_f2), index_f2=length(chdata(1).faxis);end
f=chdata(1).faxis;
temp_angle=(param.samples_per_ui*param.sample_dt)*pi.*chdata(1).faxis;
if(f(1)==0)
    temp_angle(1)=1e-20;% we don't want to divide by zero
end
PWF_tx=ones(1,length(f));
if max(upsampled_txffe) > 0
    PWF_tx=zeros(1,length(f));
    [mcur,icur] = max(upsampled_txffe);
    if exist('phase_memory','var') && ~isempty(phase_memory)
        pre_calc=1;
    else
        pre_calc=0;
    end
    for ii=1:length(upsampled_txffe)
        if upsampled_txffe(ii)==0
            %speed up:  skip cases when txffe=0
            continue;
        end
        %         PWF_tx=upsampled_txffe(ii).*exp(-1j*2*pi*(ii-icur).*f/param.fb)+PWF_tx;
        % Adee Ran 2020-06-03 remove create large 2D matrix when 1D is all
        % that is needed
        if ii==icur
            %speed up:  ii-icur=0, so just scalar addition and avoid exp calc
            PWF_tx = PWF_tx + upsampled_txffe(ii);
        else
            if pre_calc
                %speed up:  avoid vector exp calculation by externally pre-calculating it
                term_ii = upsampled_txffe(ii).*phase_memory(:,ii);
            else
                term_ii = upsampled_txffe(ii).*exp(-1j*2*pi*(ii-icur).*f/param.fb);
            end
            %bug fix:  use transpose instead of ' to avoid taking complex conjugate
            PWF_tx = PWF_tx + transpose(term_ii(:));
        end
        % /Adee
    end
end
PWF_rx=ones(1,length(f));
if exist('C','var')
    PWF_rx=zeros(1,length(f));
    for ii=-param.RxFFE_cmx:param.RxFFE_cpx
        if C(ii+param.RxFFE_cmx+1)==0
            %speed up:  skip cases when rxffe=0
            continue;
        end
        if ii+1==0
            %speed up:  ii+1=0, so just scalar addition and avoid exp calc
            PWF_rx = PWF_rx + C(ii+param.RxFFE_cmx+1);
        else
            if pre_calc
                %speed up:  avoid vector exp calculation by externally pre-calculating it
                %The latter columns of phase_memory hold RXFFE shift vectors
                term_ii=C(ii+param.RxFFE_cmx+1).*phase_memory(:,ii+param.RxFFE_cmx+1+length(upsampled_txffe));
                term_ii=transpose(term_ii);
            else
                term_ii=C(ii+param.RxFFE_cmx+1).*exp(-1j*2*pi*(ii+1).*f/param.fb);
            end
            PWF_rx=PWF_rx+term_ii;
        end
        %PWF_rx=C(ii+param.RxFFE_cmx+1).*exp(-1j*2*pi*(ii+1).*f/param.fb)+PWF_rx;
    end
end
MDFEXT=0;MDNEXT=0;MDNEXT_ICN=0;MDFEXT_ICN=0;
for ii=2:size(chdata,2)
    SINC = sin(temp_angle)./temp_angle;
    PWF_data=SINC.^2;
    PWF=PWF_data.*PWF_rx; % power weight function
    PWFnext=abs(PWF);
    PWF=PWF_data.*PWF_rx.*PWF_tx; % power weight function
    PWFfext=abs(PWF);
    if isequal(chdata(ii).type, 'FEXT')
        MDFEXT=sqrt(abs(chdata(ii).sdd21ctf).^2+MDFEXT.^2); % power sum xtk
        MDFEXT_ICN=sqrt(2*chdata(ii).delta_f/param.f2*sum( chdata(ii).A^2*PWFfext(1:index_f2).*abs(MDFEXT(1:index_f2)).^2)); %eq 46
    elseif isequal(chdata(ii).type, 'NEXT')
        MDNEXT=sqrt(abs(chdata(ii).sdd21ctf).^2+MDNEXT.^2); % power sum xtk
        MDNEXT_ICN=sqrt(2*chdata(ii).delta_f/param.f2*sum( chdata(ii).A^2*PWFnext(1:index_f2).*abs(MDNEXT(1:index_f2)).^2)); %eq 47
    end
end
if nargout == 1 && isequal(type,'NEXT')
    sigma_XT = MDNEXT_ICN*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
elseif nargout == 1 && isequal(type,'FEXT')
    sigma_XT = MDFEXT_ICN*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
elseif nargout == 3
    sigma_XT=norm([ MDNEXT_ICN MDFEXT_ICN ])*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
    sigma_NEXT = MDNEXT_ICN*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
    sigma_FEXT = MDFEXT_ICN*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
end

function out=hrem(h,index,N_bf,bmaxg)

out=[ h(1:index-1) h(index:index+N_bf-1)- sign(h(index:index+N_bf-1)).*   (min( bmaxg, abs( h(index:index+N_bf-1) )))  h(index+N_bf:end) ];
% faster than single line function
% hrem =@(h,index,N_bf,bmaxg) ...
%     [ h(1:index-1) ...
%     h(index:index+N_bf-1)- sign(h(index:index+N_bf-1)).*   (min( bmaxg, abs( h(index:index+N_bf-1) ))) ...
%     h(index+N_bf:end) ]...
%     ;

%% floating DFE taps
function [Sout] = interp_Sparam(Sin,fin,fout, ...
    opt_interp_Sparam_mag, opt_interp_Sparam_phase,OP)
% Sout = interp_Sparam(Sin,fin,fout)
%
% Interpolate S-parameters Sin from frequency grid fin to frequency grid
% fout.

if ( fin(end)<fout(end) )
    %    warning('Channel high frequencies extrapolation might be inaccurate!');
end

H_mag = abs(Sin);
H_mag(H_mag<eps)=eps; % handle ill cases...
H_ph = unwrap(angle(Sin));
% For long delay channels, the result can turn anti-causal if frequency step is too coarse. Don't let the
% user ignore that.
if mean(diff(H_ph))>0
    if OP.DEBUG
        warning('Anti-causal response found. Finer frequency step is required for this channel');
    else
        error('Anti-causal response found. Finer frequency step is required for this channel');
    end
end

%opt_interp_Sparam_mag='linear_trend_to_DC';
switch opt_interp_Sparam_mag
    case {'linear_trend_to_DC','linear_trend_to_DC_log_trend_to_inf'}
        if -iscolumn(H_mag), H_mag=H_mag.';end
        if -iscolumn(fin), fin=fin.';end
        fin_x=fin;
        H_mag_x=H_mag(:);
        if fin(1)>0
            p=polyfit(fin(1:10), H_mag(1:10), 1);
            dc_trend_val=polyval(p, 0);
            fin_x=[0, fin_x];
            H_mag_x = [dc_trend_val; H_mag_x];
        end
        if fin(end)<fout(end)
            warn_state=warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
            mid_freq_ind=round(length(fin)/2);
            p=polyfit(fin(mid_freq_ind:end), H_mag(mid_freq_ind:end), 1);
            warning(warn_state);
            hf_trend_val=polyval(p, fout(end));
            if hf_trend_val>H_mag(end)
                hf_trend_val=H_mag(end);
                hf_logtrend_val = H_mag(end);
            elseif hf_trend_val<eps
                hf_trend_val=eps;
                hf_logtrend_val = realmin;
            end
            fin_x=[fin_x, fout(end)];
            H_mag_x = [H_mag_x; hf_trend_val];
        end
        H_mag_i = interp1(fin_x, H_mag_x, fout, 'linear', 'extrap');
        if strcmp(opt_interp_Sparam_mag,'linear_trend_to_DC_log_trend_to_inf')
            H_logmag_i = exp(interp1(fin_x, log([H_mag_x(1:end-1);hf_logtrend_val]), fout, 'linear', 'extrap'));
            indx = find(fout > fin(end),1,'first');
            H_mag_i(indx:end) = H_logmag_i(indx:end);
        end
    case 'trend_to_DC'
        % extrapolate to trend value at DC.
        if -iscolumn(H_mag), H_mag=H_mag.';end
        if -iscolumn(fin), fin=fin.';end
        fin_x=fin;
        H_mag_x=H_mag;
        if fin(1)>0
            warn_state=warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
            p=polyfit(fin(1:10), log10(H_mag(1:10)), 1);
            dc_trend_val=10^polyval(p, 0);
            fin_x=[0, fin_x];
            H_mag_x = [dc_trend_val H_mag_x];
        end
        if fin(end)<fout(end)
            warn_state=warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
            mid_freq_ind=round(length(fin)/2);
            p=polyfit(fin(mid_freq_ind:end), log10(H_mag(mid_freq_ind:end)), 1);
            warning(warn_state);
            hf_trend_val=10^polyval(p, fout(end));
            if hf_trend_val>H_mag(end)
                hf_trend_val=H_mag(end);
            end
            fin_x=[fin_x, fout(end)];
            H_mag_x = [H_mag_x hf_trend_val];
        end
        H_mag_i = 10.^interp1(fin_x,log10(H_mag_x),fout,'linear', 'extrap');
    case 'extrap_to_DC_or_zero'
        % same as extrap_to_DC but detect AC-coupled channels and
        % extrapolate them to 0.
        if fin(1)>0 && 20*log10(H_mag(1))<-20
            % assume AC coupling, with 0 at DC
            H_mag_i = 10.^interp1([0, fin],[-100; log10(H_mag)],fout(fout<=fin(end)),'linear', 'extrap');
        else
            H_mag_i = 10.^interp1(fin, log10(H_mag), fout(fout<=fin(end)),'linear', 'extrap');
        end
        H_mag_i(fout>fin(end)) = H_mag(end);
    case 'extrap_to_DC'
        % first extrapolate down to DC, then use highest available frequency
        % for higher frequencies
        H_mag_i = 10.^interp1(fin,log10(H_mag),fout(fout<=fin(end)),'linear', 'extrap');
        H_mag_i(fout>fin(end)) = H_mag(end);
    case 'old'
        H_mag_i = interp1(fin,H_mag,fout,'linear','extrap');
    otherwise
        error('COM:Extrap:InvalidOption', 'opt_interp_Sparam_mag valid values are "old", "extrap_to_DC"');
end

H_ph_i = interp1(fin,squeeze(H_ph),fout,'linear', 'extrap');

%opt_interp_Sparam_phase='trend_and_shift_to_DC';
%opt_interp_Sparam_phase='interp_cubic_to_dc_linear_to_inf';
switch opt_interp_Sparam_phase
    case 'old'
        H_ph_i = H_ph_i-H_ph_i(1);
    case 'zero_DC'
        H_ph_i(1) = 0;
    case 'interp_to_DC'
        if fin(1) ~= 0
            H_ph_i = interp1([0; fin(:)], [0; H_ph(:)], fout, 'linear', 'extrap');
        end
    case 'extrap_cubic_to_dc_linear_to_inf'
        if fin(1) ~= 0
            % estimate low frequency group delay
            group_delay = -diff(H_ph(:))./diff(fin(:));
            low_freq_gd = group_delay(1:50);
            %  calculate trend, throwing away outliers
            m = median(low_freq_gd); sigma = std(low_freq_gd);
            lf_trend = mean(low_freq_gd(abs(low_freq_gd-m)<sigma));
            % correct outliers in first 10 phase samples
            for k=10:-1:1
                H_ph(k) = H_ph(k+1) + lf_trend*(fin(k+1)-fin(k));
            end
            H_ph_cubic = interp1(fin, H_ph, fout, 'pchip', 'extrap');
            H_ph_linear = interp1(fin, H_ph, fout, 'linear', 'extrap');
            % modification - trend to inf
            if (1)
                high_freq_gd = group_delay(end-50:end);
                %  calculate trend, throwing away outliers
                m = median(high_freq_gd); sigma = std(high_freq_gd);
                hf_trend = -mean(high_freq_gd(abs(high_freq_gd-m)<sigma));
                hf_extrap_range = find(fout>fin(end));
                last_data_sample = hf_extrap_range(1)-1;
                H_ph_linear(hf_extrap_range) = H_ph_linear(last_data_sample) + (fout(hf_extrap_range)-fout(last_data_sample))*hf_trend;
                %                for k=hf_range
                %                    H_ph_linear(k) = H_ph_linear(k-1) + hf_trend*(fout(k)-fout(k-1));
                %                end
            end
            
            [UNUSED_OUTPUT, indx] = min(abs(H_ph_cubic-H_ph_linear)); %#ok<ASGLU>
            H_ph_i=H_ph_cubic;
            H_ph_i(indx:end) = H_ph_linear(indx:end);
            H_ph_i = H_ph_linear; % John Ewen 12/13/2019
        end
    case 'interp_and_shift_to_DC'
        if fin(1) ~= 0
            dc_phase_trend = H_ph(1)-(H_ph(2)-H_ph(1))/(fin(2)-fin(1))*fin(1);
            H_ph_i = interp1([0; fin(:)], [0; H_ph(:)-dc_phase_trend], fout, 'linear', 'extrap');
        end
    case 'trend_and_shift_to_DC'
        % estimate low frequency group delay
        group_delay = -diff(H_ph(:))./diff(fin(:));
        low_freq_gd = group_delay(1:50);
        %  calculate trend, throwing away outliers
        m = median(low_freq_gd); sigma = std(low_freq_gd);
        lf_trend = mean(low_freq_gd(abs(low_freq_gd-m)<sigma));
        fin_x=fin;
        H_ph_x=H_ph(:);
        if fin(1) ~= 0
            % correct outliers in first 10 phase samples
            for k=10:-1:1
                H_ph(k) = H_ph(k+1) + lf_trend*(fin(k+1)-fin(k));
            end
            
            % shift all phase data so that DC extrapolation to 0 follows trend
            dc_phase_trend = H_ph(1)+lf_trend*(fin(1)-0);
            fin_x=[0, fin_x];
            H_ph_x=[0; H_ph(:)-dc_phase_trend];
        end
        % Modification: extrapolate using trend. (interp1 with "extrap" extrapolates using just
        % the last two samples, so noise can create an inverted slope and
        % non-causal response).
        if fout(end)>fin(end)
            group_delay = -diff(H_ph_x(:))./diff(fin_x(:));
            % p=polyfit(fin_x', H_ph_x, 1);
            hf_phase_trend = H_ph_x(end)-median(group_delay)*(max(fout)-max(fin_x));
            % hf_phase_trend=polyval(p,max(fout));
            fin_x=[fin_x, fout(end)];
            H_ph_x=[H_ph_x; hf_phase_trend];
        end
        H_ph_i = interp1(fin_x, H_ph_x, fout, 'linear', 'extrap');
        
    otherwise
        error('COM:Extrap:InvalidOption', ...
            'debug_interp_Sparam valid values are "old", "zero_DC", "interp_to_DC", "interp_and_shift_to_DC", "trend_and_shift_to_DC", "interp_cubic_to_dc_linear_to_inf"');
end
H_i = H_mag_i.*exp(1j*H_ph_i);
Sout=H_i;
function [ s11out, s12out, s21out, s22out]=make_full_pkg(type,faxis,param,channel_type,mode,include_die)

%This function makes the TX or RX package.  The type input must be
%'TX' or 'RX'
%If the mode argument is omitted, mode='dd' is assumed.  Currently
%mode='dc' is only used when making the TX package for AC CM noise
%inclusion.  The Rx package for 'dc' mode is still generated using
%the same parameters as 'dd' mode
%channel_type should be 'THRU' 'FEXT' or 'NEXT'
%
%One instance of package block looks like this (if no elements are set to 0):
%-------------Lcomp----------Tline---------------
%   |                   |               |
%   Cpad                Cbump           Cball
%   |                   |               |
%------------------------------------------------

if nargin<6
    %optional input "include_die"=0 allows die parameters to be forced to 0
    %this includes Cpad, Lcomp, and Cbump
    include_die=1;
end
if nargin<5
    mode='dd';
end


if ~isempty(param.PKG_NAME)
    %The gamma and tau parameters do not currently have a separate Tx and Rx home to live (they were locked for both sides originally)
    %so they are swapped in depending on if Tx or Rx is set for type
    %Note that param is not returned from this function, so the swap does not persist
    swap_fields = {'pkg_gamma0_a1_a2' 'pkg_tau'};
    if strcmpi(type,'tx')
        pkg_name = param.PKG_NAME{1};
    elseif strcmpi(type,'rx')
        pkg_name = param.PKG_NAME{2};
    else
        error('Pkg type must be Tx or Rx');
    end
    pkg_parameter_struct = param.PKG.(pkg_name);

    
    for j=1:length(swap_fields)
        param.(swap_fields{j}) = pkg_parameter_struct.(swap_fields{j});
    end
    
end

C_diepad = param.C_diepad;
C_pkg_board = param.C_pkg_board;
% [ahealey] Unpack optional compensating L and "bump" C model parameters.
L_comp = param.L_comp;
C_bump = param.C_bump;
if ~include_die
    %best to multiply by 0.  that way vectors maintain original size
    C_diepad=C_diepad*0;
    L_comp=L_comp*0;
    C_bump=C_bump*0;
end
% [ahealey] End of modifications.
% generate TX package according to channel type.
[ncases, mele]=size(param.z_p_next_cases);

%Syntax update for C_diepad and L_comp
%Allow a chain of values to be entered as a matrix:
%[L_Tx1 L_Tx2 L_Tx3 ; L_Rx1 L_Rx2 L_Rx3]
if isvector(C_diepad)
    Cd_Tx=C_diepad(1);
    Cd_Rx=C_diepad(2);
    L_comp_Tx=L_comp(1);
    L_comp_Rx=L_comp(2);
    num_blocks=mele;
else
    Cd_Tx=C_diepad(1,:);
    Cd_Rx=C_diepad(2,:);
    L_comp_Tx=L_comp(1,:);
    L_comp_Rx=L_comp(2,:);
    num_blocks=mele+length(Cd_Tx)-1;
end
extra_LC=length(Cd_Tx)-1;
%note:  "insert_zeros" is empty if length(Cd_Tx) = 1
insert_zeros=zeros([1 extra_LC]);

%Updated technique of building Tx/Rx packages
%each index corresponds to the package segment
switch type
    case 'TX'
        switch mele
            case 1
                Cpad=Cd_Tx;
                Lcomp=L_comp_Tx;
                Cbump=C_bump(1);
                Cball=C_pkg_board(1);
                Zpkg=param.pkg_Z_c(1);
            case 4
                Cpad=[Cd_Tx 0 0 0];
                Lcomp=[L_comp_Tx 0 0 0];
                Cbump=[C_bump(1) 0 0 0];
                Cball=[0 0 param.C_v(1) C_pkg_board(1)];
                Zpkg=param.pkg_Z_c(1,:);
            otherwise
                error('package syntax error')
        end
        switch upper(channel_type)
            case 'THRU'
                Len=param.Pkg_len_TX;
            case 'NEXT'
                Len=param.Pkg_len_NEXT;
            case 'FEXT'
                Len=param.Pkg_len_FEXT;
        end
    case 'RX'
        switch mele
            case 1
                Cpad=Cd_Rx;
                Lcomp=L_comp_Rx;
                Cbump=C_bump(2);
                Cball=C_pkg_board(2);
                Zpkg=param.pkg_Z_c(2);
            case 4
                Cpad=[Cd_Rx 0 0 0];
                Lcomp=[L_comp_Rx 0 0 0];
                Cbump=[C_bump(2) 0 0 0];
                Cball=[0 0 param.C_v(2) C_pkg_board(2)];
                Zpkg=param.pkg_Z_c(2,:);
            otherwise
                error('package syntax error')
        end
        switch upper(channel_type)
            case 'THRU'
                Len=param.Pkg_len_RX;
            case 'NEXT'
                Len=param.Pkg_len_RX;
            case 'FEXT'
                Len=param.Pkg_len_RX;
        end
end

%Insert the extra 0 at the front end of Cball, Cbump, Len, and Zpkg
Cball=[insert_zeros Cball];
Cbump=[insert_zeros Cbump];
Len=[insert_zeros Len];
Zpkg=[insert_zeros Zpkg];

% debug_string='';
% for j=1:length(Zpkg)
%     if Cpad(j)~=0
%         debug_string=[debug_string sprintf(', Cd=%0.4g',Cpad(j))];
%     end
%     if Lcomp(j)~=0
%         debug_string=[debug_string sprintf(', Ls=%0.4g',Lcomp(j))];
%     end
%     if Cbump(j)~=0
%         debug_string=[debug_string sprintf(', Cb=%0.4g',Cbump(j))];
%     end
%     if Len(j)~=0
%         debug_string=[debug_string sprintf(', Len=%0.4g Zc=%0.3g',Len(j),Zpkg(j))];
%     end
%     if Cball(j)~=0
%         debug_string=[debug_string sprintf(', Cp=%0.4g',Cball(j))];
%     end
% end
% if length(debug_string)>2
%     debug_string=debug_string(3:end);
% end

% tx package
pkg_param=param;
if strcmpi(mode,'dc')
    % change tx package to CC mode
    pkg_param.Z0=pkg_param.Z0/2;
    Cpad=Cpad*2;
    Cball=Cball*2;
    Zpkg=Zpkg*2;
    Lcomp=Lcomp/2;
    Cbump=Cbump*2;
end
switch num_blocks
    case 1
        [ s11out, s12out, s21out, s22out ]= make_pkg(faxis, Len(1), Cpad(1), Cball(1),Zpkg(1), pkg_param, Lcomp(1), Cbump(1));
    otherwise
        for j=1:num_blocks
            [spkg11,spkg12,spkg21,spkg22]=make_pkg(faxis, Len(j),  Cpad(j),Cball(j) ,Zpkg(j), pkg_param, Lcomp(j),Cbump(j));
            if j==1
                s11out=spkg11; s12out=spkg12; s21out=spkg21; s22out=spkg22;
            else
                [ s11out, s12out, s21out, s22out ]=combines4p(  s11out, s12out, s21out, s22out, spkg11,spkg12,spkg21,spkg22   );
            end
        end
end
function [ s11out, s12out, s21out, s22out ] = make_pkg(f, pkg_len, cpad, cball, pkg_z, param, varargin)
f(f<eps)=eps;
tau = param.pkg_tau; gamma0_a1_a2=param.pkg_gamma0_a1_a2; Lenscale=pkg_len; zref=param.Z0;
%% Equation 93A-8
s11pad= -1i*2*pi.*f*cpad*zref./(2+1i*2*pi.*f*cpad*zref);
s21pad= 2./(2+1i*2*pi.*f*cpad*zref);

% [ahealey] Add compensating L and shunt C (bump) when requested.
s12pad = s21pad;
s22pad = s11pad;
if nargin > 6
    lcomp = varargin{1};
    if lcomp>0
        s11comp = (1i*2*pi*f*lcomp/zref)./(2+1i*2*pi*f*lcomp/zref);
        s21comp = 2./(2+1i*2*pi*f*lcomp/zref);
        [s11pad, s12pad, s21pad, s22pad] = combines4p( ...
            s11pad, s12pad, s21pad, s22pad, ...
            s11comp, s21comp, s21comp, s11comp);
    end
end
if nargin > 7
    cbump = varargin{2};
    if cbump>0
        s11bump = -1i*2*pi.*f*cbump*zref./(2+1i*2*pi.*f*cbump*zref);
        s21bump = 2./(2+1i*2*pi.*f*cbump*zref);
        [s11pad, s12pad, s21pad, s22pad] = combines4p( ...
            s11pad, s12pad, s21pad, s22pad, ...
            s11bump, s21bump, s21bump, s11bump);
    end
end
% [ahealey] End of modifications.

[ S11, S12, S21, S22 ] = synth_tline(f, pkg_z, zref, gamma0_a1_a2, tau, Lenscale); %#ok<NASGU,ASGLU>
% [ahealey] Symmetry cannot be assumed with more complex termination models.
% [ s11out1, s12out1, s21out1, s22out1 ]= ...
%     combines4p(  s11pad, s21pad, s21pad, s11pad, S11, S21, S21, S11 ); % first part of equation 93A-15
[s11out1, s12out1, s21out1, s22out1] = combines4p( ...
    s11pad, s12pad, s21pad, s22pad, ...
    S11, S21, S21, S11);
% [ahealey] End of modifications.

%% Equation 93A-8
s11ball= -1i*2*pi.*f*cball*zref./(2+1i*2*pi.*f*cball*zref);
s21ball= 2./(2+1i*2*pi.*f*cball*zref);
[ s11out, s12out, s21out, s22out ]= ...
    combines4p( s11out1, s12out1, s21out1, s22out1, s11ball, s21ball, s21ball, s11ball );% second part of equation 93A-15

function missingParameter (parameterName)
error( 'error:badParameterInformation', ...
'The data for mandatory parameter %s is missing or incorrect' , parameterName);

function pdf = normal_dist(sigma,nsigma,binsize)
pdf.BinSize=binsize;
pdf.Min=-round(2*nsigma*sigma/binsize); % RIM 03/03/2023 capture more of the tails
pdf.x=(pdf.Min:-pdf.Min)*binsize;
pdf.y=exp(-pdf.x.^2/(2*sigma^2+eps));
pdf.y=pdf.y/sum(pdf.y);

function result=optimize_fom(OP, param, chdata, sigma_bn,do_C2M)
%% input
% chdata(1).uneq_imp_response is the impulse response input expected to be normalized to At, peak drive voltage
% baud_rate - baud rate in seconds
% param.samples_per_ui = samples per UI of IR
% param.max_ctle - maximum ac to dc gain in dB
% param.tx_ffe(1) - maximum pre cursor (positive value)
% param.tx_ffe(2) - maximum post cursor (positive value)
% param.tx_ffe_step - sweep step size for tx pre and post taps
% param.ndfe - number of reference dfe taps
% do_C2M.  set to 0 for standard optimize_fom.  set to 1 for optimize_fom_for_C2M
% output
% result.eq.txle - [ precusor curosr postcursor]: pre and post are negative
% result.eq.ctle - index of CTLE parameters in table
% result.IR - impulse response
% result.avail_signal - maximum signal after equalization
% result.avail_sig_index - index in result.IR of max signal
% result.best_FOM - best raw ISI


min_number_of_UI_in_response=40;
baud_rate=1/param.ui;
% H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(1).faxis/(param.f_r*param.fb));
f=chdata(1).faxis;

%Read user input of ts_sample_adj_range
%if one value was entered, go from 0 to that value
%if 2 values were entered, go from the 1st value to the 2nd value
if length(param.ts_sample_adj_range)==1
    param.ts_sample_adj_range(2)=param.ts_sample_adj_range(1);
    param.ts_sample_adj_range(1)=0;
end
full_sample_range=param.ts_sample_adj_range(1):param.ts_sample_adj_range(2);

H_bt=Bessel_Thomson_Filter(param,f,OP.Bessel_Thomson);
H_bw=Butterworth_Filter(param,f,OP.Butterworth);
H_RCos=Raised_Cosine_Filter(param,f,OP.Raised_Cosine);% experiment with RCos
% need to include H_RCos in noise and when computing the system ir for thru
% and crosstalk
H_r=H_bw.*H_bt.*H_RCos;
%% Bill Kirkland, need to get auto correlation of H_r.*HCTLE
% Get f vector from 0 to Fs/2-delta_f.
N_fft_by2 = 512;
f_xc = ((1:N_fft_by2)-1)/N_fft_by2*baud_rate*param.samples_per_ui/2;
H_bt_xc=Bessel_Thomson_Filter(param,f_xc,OP.Bessel_Thomson);
H_bw_xc=Butterworth_Filter(param,f_xc,OP.Butterworth);
H_RCos_xc=Raised_Cosine_Filter(param,f_xc,OP.Raised_Cosine);
H_r_xc=H_bw_xc.*H_bt_xc.*H_RCos_xc;
%%

% system noise H_sy PSD
if OP.USE_ETA0_PSD
    fspike=1e9;
    % requires communication tool box if used
    H_sy=sinc(sqrt(2)*(chdata(1).faxis-fspike)/(fspike)).^2;
else
    H_sy=ones(1,length(chdata(1).faxis));
end

%Build txffe values dynamically
%any param field that is "tx_ffe_cm<X>_values" is a precursor
%any param field that is "tx_ffe_cp<X>_values" is a postcursor
%where <X> is any integer
param_fields=fieldnames(param);
num_pre=length(find(~cellfun('isempty',regexp(param_fields,'tx_ffe_cm\d+_values'))));
num_post=length(find(~cellfun('isempty',regexp(param_fields,'tx_ffe_cp\d+_values'))));
num_taps=num_pre+num_post;
cur=num_pre+1;
%txffe_cell combines all the txffe values into a single cell array
%It is ordered: [precursorN precursorN-1 ... precursor1 postcursor1 ... postcursorN-1 postcursorN]
txffe_cell=cell(1,num_taps);
for k=num_pre:-1:1
    idx=num_pre-k+1;
    this_tx_field=sprintf('tx_ffe_cm%d_values',k);
    txffe_cell{idx}=param.(this_tx_field);
end
for k=1:num_post
    idx=k+num_pre;
    this_tx_field=sprintf('tx_ffe_cp%d_values',k);
    txffe_cell{idx}=param.(this_tx_field);
end
%total number of txffe runs is the product of the lengths of each tap
txffe_lengths=cellfun('length',txffe_cell);
if isempty(txffe_cell)
    num_txffe_runs=1;
else
    num_txffe_runs=prod(txffe_lengths);
end
%txffe_sweep_indices are used in the LOCAL_SEARCH block
%any tap with length=1 can be ignored
%Also is statistically likely that taps with greater number of values
%will exceed the LOCAL SEARCH criteria, so searching those first is faster
txffe_sweep_indices=find(txffe_lengths>1);
[~,length_sort]=sort(txffe_lengths(txffe_sweep_indices),'descend');
txffe_sweep_indices=txffe_sweep_indices(length_sort);
num_txffe_sweep_indices=length(txffe_sweep_indices);

gdc_values = param.ctle_gdc_values;
Gffe_values = param.cursor_gain;
switch param.CTLE_type
    case 'CL93'
    case 'CL120d'
        g_DC_HP_values =param.g_DC_HP_values;
    case 'CL120e'
        f_HP_Z=param.f_HP_Z;
        f_HP_P=param.f_HP_P;
        
end
best_ctle = [];
best_FOM = -inf;
best_txffe = [];
delta_sbr = [];
PSD_results=[];
MMSE_results=[];
best_bmax=param.bmax;
%AJG021820
best_bmin=param.bmin;
h_J=[];
pxi=0;
if OP.DISPLAY_WINDOW
    hwaitbar=waitbar(0);
else
    fprintf('FOM search ');
end
FOM=0;
if ~OP.RxFFE
    Gffe_values=0;
end
param.ndfe_passed=param.ndfe;
old_loops=0;
new_loops=0;

%GDC Qual construction
gqual= param.gqual;
g2qual=param.g2qual;
if ~strcmp(param.CTLE_type,'CL120d')
    qual=ones(1,length(gdc_values));
else
    if isempty(gqual) && isempty(g2qual)
        qual=ones(length(g_DC_HP_values),length(gdc_values));
    else
        qual=zeros(length(g_DC_HP_values),length(gdc_values));
        
        %prepare gqual and g2qual
        [g2qual,si]=sort(g2qual,'descend');
        gqual=gqual(si,:);
        tmp=g2qual;
        g2qual=zeros(length(tmp),2);
        for kk=1:length(tmp)
            if kk==1
                g2qual(kk,:)=[tmp(kk)+eps tmp(kk)];
            else
                g2qual(kk,:)=[tmp(kk-1) tmp(kk)];
            end
            gqual(kk,:)=sort(gqual(kk,:),'descend');
        end
        
        %Qual Construction
        for jj=1:length(g_DC_HP_values)
            for ii=1:length(gdc_values)
                for kk=1:size(gqual,1)
                    if g_DC_HP_values(jj) >= g2qual(kk,2) && g_DC_HP_values(jj) < g2qual(kk,1)
                        if gdc_values(ii) >= gqual(kk,2) && gdc_values(ii) < gqual(kk,1)
                            qual(jj,ii)=1;
                            break;
                        end
                    end
                end
            end
        end
    end
end

progress_interval=0.025;
if do_C2M
    loop_count=[1 2];
    T_O=floor((param.T_O/1000)*param.samples_per_ui);
    T_O=max(0,T_O);
else
    loop_count=1;
    T_O=0;
end
switch param.CTLE_type
    case 'CL93'
        lf_indx=1;
    case 'CL120d'
        lf_indx=length(g_DC_HP_values);
    case 'CL120e'
        lf_indx=1;
end
runs=length(gdc_values)*lf_indx*length(Gffe_values)*num_txffe_runs;
if  OP.Optimize_loop_speed_up == 1
    OP.BinSize = 1e-4;
    OP.impulse_response_truncation_threshold = 1e-3;
end

%Used to speed up FFE by only performing circshift when necessary
pulse_ctle_circshift=[];
ctle_response_updated=1;

%Used to speed up get_xtlk_noise by pre-calculating all the phase shift exponentials
calc_exp_phase=0;

%calculate cur index and pre/post indices outside of the loop
cur_start=cur;
precursor_indices=[];
postcursor_indices=[];
auto_count_trigger=0;
for kv=1:num_taps
    if ~auto_count_trigger && length(txffe_cell{kv})==1 && txffe_cell{kv}==0
        %precursor values fill the beginning of the vector.  Any empty precursor means
        %cursor position must be subtracted by 1
        if kv<cur_start
            cur=cur-1;
        end
    else
        %non empty value:  add to precursor or postcursor indices depending on position
        %in the vector
        if kv<cur_start
            auto_count_trigger=1;
            precursor_indices=[precursor_indices kv];
        else
            auto_count_trigger=0;
            postcursor_indices=[postcursor_indices kv];
        end
    end
end
if ~isempty(postcursor_indices)
    postcursor_indices=postcursor_indices(1):postcursor_indices(end);
end

%Calculate the full grid matrix of all txffe combinations
if isempty(txffe_cell)
    TXFFE_grid=0;
    FULL_tx_index_vector=1;
else
    TXFFE_grid=Full_Grid_Matrix(txffe_cell);
    %Also calculate the full grid matrix for the index used in each txffe combination
    %(the index is used in the LOCAL SEARCH block)
    for k=1:num_taps
        txffe_index_cell{k}=1:txffe_lengths(k);
    end
    FULL_tx_index_vector=Full_Grid_Matrix(txffe_index_cell);
end

%pre-calculate cursor to save time
txffe_cursor_vector=1-sum(abs(TXFFE_grid),2);

%pre-calculate full txffe for each iteration to save time
precursor_matrix=TXFFE_grid(:,precursor_indices);
postcursor_matrix=TXFFE_grid(:,postcursor_indices);
txffe_matrix = [precursor_matrix txffe_cursor_vector  postcursor_matrix];

if OP.TDMODE
    uneq_field='uneq_pulse_response';
    ctle_field='ctle_pulse_response';
else
    uneq_field='uneq_imp_response';
    ctle_field='ctle_imp_response';
end

%Speed up search for max(sbr)
if OP.TDMODE
    [~,init_max]=max(chdata(1).uneq_pulse_response);
else
    [~,init_max]=max(filter(ones(1,param.samples_per_ui),1,chdata(1).uneq_imp_response));
end
UI_max_window=20;
start_max_idx=init_max-UI_max_window*param.samples_per_ui;
if start_max_idx<1
    start_max_idx=1;
end
end_max_idx=init_max+UI_max_window*param.samples_per_ui;
if end_max_idx>length(chdata(1).(uneq_field))
    end_max_idx=length(chdata(1).(uneq_field));
end

itick_skips=0;
itick_cases=0;
FOM_TRACKER(1:length(Gffe_values),1:length(gdc_values),1:lf_indx,1:size(TXFFE_grid,1),1:length(full_sample_range))=0;
for i=loop_count
    
    for Gffe_index=1:length(Gffe_values)
        param.current_ffegain=Gffe_values(Gffe_index);
        for ctle_index=1:length(gdc_values)
            g_dc = gdc_values(ctle_index);
            kacdc = 10^(g_dc/20);
            CTLE_fp1 = param.CTLE_fp1(ctle_index);
            CTLE_fp2 = param.CTLE_fp2(ctle_index);
            CTLE_fz = param.CTLE_fz(ctle_index);
            switch param.CTLE_type
                case 'CL93'
%
                case 'CL120d'
%
                case 'CL120e'
                    HP_Z = param.f_HP_Z(ctle_index);
                    HP_P = param.f_HP_P(ctle_index);
            end
              %% HF Boost
            ctle_gain = (kacdc + 1i*chdata(1).faxis/CTLE_fz) ./ ...
                ((1+1i*chdata(1).faxis/CTLE_fp1).*(1+1i*chdata(1).faxis/CTLE_fp2));
%% Mid Frequency Boost
            ctle_gain_xc = (kacdc + 1i*f_xc/CTLE_fz) ./ ...
                ((1+1i*f_xc/CTLE_fp1).*(1+1i*f_xc/CTLE_fp2)); % Bill Kirkland          
            for  g_LP_index=1:lf_indx
                
                %GDC Qual Check
                if qual(g_LP_index,ctle_index)==0
                    pxi=pxi+num_txffe_runs;
                    continue;
                end
                
                switch param.CTLE_type
                    case 'CL93'
                        H_low=1;
                        kacde_DC_low=1;
                    case 'CL120d'
                        g_DC_low = g_DC_HP_values(g_LP_index);
                        f_HP=param.f_HP(g_LP_index);
                        kacde_DC_low = 10^(g_DC_low/20);
                        H_low=(kacde_DC_low +  1i*chdata(1).faxis/f_HP)./(1 + 1i*chdata(1).faxis/f_HP);
                        H_low_xc = (kacde_DC_low +  1i*f_xc/f_HP)./(1 + 1i*f_xc/f_HP);% Bill Kirkland
                    case 'CL120e' % z1 has been adusted on read in
                        H_low=(1 +  1i*chdata(1).faxis/HP_Z)./(1 + 1i*chdata(1).faxis/HP_P);
                        H_low_xc=(1 +  1i*f_xc/HP_Z)./(1 + 1i*f_xc/HP_P); % Bill Kirkland
                end
                H_ctf=H_low.*ctle_gain;
                switch upper(OP.FFE_OPT_METHOD)
                    case 'WIENER-HOPF'
                                        %% Bill Kirkland
                        H_ctf_xc = H_low_xc.*ctle_gain_xc;
                        H_rx_ctle_xc = H_r_xc.*H_ctf_xc;
                        % use Fourier Transform pair for correlation as we have to
                        % take ifft of H_r anyways.
                        % onesided and two sided responses - tricky, tricky, tricky
                        Var_eta0 =  param.eta_0*f_xc(end)/1e9;
                        XC_rx_ctle = ifft (H_rx_ctle_xc.*conj(H_rx_ctle_xc),2*length(H_rx_ctle_xc),'symmetric');
                        Noise_XC = Var_eta0.*XC_rx_ctle(1:param.samples_per_ui:N_fft_by2);

                        if OP.Do_White_Noise
                            Noise_XC = Noise_XC(1);
                        end
                    otherwise
                        Noise_XC=[];            
                end


                
                if OP.INCLUDE_CTLE==1
                    for k=1:param.num_s4p_files
                        ir_peak = max(abs(chdata(k).(uneq_field)));
                        ir_last  = find(abs(chdata(k).(uneq_field))>ir_peak*OP.impulse_response_truncation_threshold, 1, 'last');
                        chdata(k).(uneq_field) = chdata(k).(uneq_field)(1:ir_last);
                        chdata(k).(ctle_field) = TD_CTLE(chdata(k).(uneq_field), baud_rate ...
                            , CTLE_fz, CTLE_fp1, CTLE_fp2, g_dc, param.samples_per_ui);
                        switch param.CTLE_type
                            case 'CL93'
                            case 'CL120d'
                                chdata(k).(ctle_field) = TD_CTLE(chdata(k).(ctle_field), baud_rate, f_HP, f_HP,100e100 , g_DC_low , param.samples_per_ui);
                            case 'CL120e' % z1 has been adusted on read in
                                chdata(k).(ctle_field) = TD_CTLE(chdata(k).(ctle_field), baud_rate, HP_Z,HP_P,100e100 , 0 , param.samples_per_ui);
                        end
                    end
                    %set the flag to show ctle response was updated
                    ctle_response_updated=1;
                else
                    for k=1:param.num_s4p_files
                        chdata(k).(ctle_field) = chdata(k).(uneq_field);
                    end
                end
                for k=1:param.num_s4p_files
                    chdata(k).sdd21ctf=chdata(k).sdd21.*H_ctf; % sdd21 is a VTF, includes H_t, H_f, and package
                end
                %% Equation 93A-22 %%
                %         figure(1000)
                %         semilogx(chdata(1).faxis/1e9,db(H_ctf))
                %         hold on
                if OP.RX_CALIBRATION
                    ctle_gain2 = (kacdc + 1i*chdata(2).faxis/CTLE_fz) ./ ...
                        ((1+1i*chdata(2).faxis/CTLE_fp1).*(1+1i*chdata(2).faxis/CTLE_fp2));
                    switch param.CTLE_type
                        case 'CL93'
                            H_low2=1;
                        case 'CL120d'
                            g_DC_low = g_DC_HP_values(g_LP_index);
                            f_HP=param.f_HP(g_LP_index);
                            kacde_DC_low = 10^(g_DC_low/20);
                            H_low2=(kacde_DC_low +  1i*chdata(1).faxis/f_HP)./(1 + 1i*chdata(1).faxis/f_HP);
                        case 'CL120e' % z1 has been adusted on read in
                            H_low2=(1 +  1i*chdata(1).faxis/HP_Z)./(1 + 1i*chdata(1).faxis/HP_P);
                    end
                    H_ctf2=H_low2.*ctle_gain2;
                end
                % RIM 11-30-2020 moved to a subfunction
                [sigma_N] = get_sigma_eta_ACCM_noise(chdata,param,H_sy,H_r,H_ctf);
                if OP.RX_CALIBRATION
                    sigma_ne = get_sigma_noise( H_ctf2, param, chdata, sigma_bn); %% Equation 93A-48 %%
                    sigma_NEXT=sqrt(param.eta_0*sum( abs(H_sy(2:end).^2 .* H_r(2:end).*2 .* H_ctf(2:end).^2 ) .* diff(chdata(1).faxis)/1e9));% changed from /chdata(1).faxis(end) B. Kirkland S. Elnagar 11/6/2021
                else
                    %% Equations 93A-33 and 93A-34  for NEXT - independent of TXFFE setting %%
                    sigma_NEXT =  get_xtlk_noise( [0 1 0], 'NEXT', param, chdata );
                    sigma_ne=0;
                end
                
                if param.GDC_MIN ~= 0 && gdc_values(ctle_index) + g_DC_HP_values(g_LP_index) > param.GDC_MIN
                    pxi=pxi+num_txffe_runs;
                    continue; % change per 0.3k draft 2.3
                end
                %%
                PSD_results=[];
                if strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE
                    OP.WO_TXFFE=1;
                    PSD_results=get_PSDs(PSD_results,[],[],[],gdc_values(ctle_index),g_DC_low,param,chdata,OP);
                end
                %TXFFE Loop
                %Originally this was a separate for loop for each tap, but it is now all contained in the TXFFE_grid matrix to use a single modular loop
                for TK=1:size(TXFFE_grid,1)

                    pxi=pxi+1;
                    progress = pxi/runs;
                    if OP.DISPLAY_WINDOW
                        if ~mod(pxi,floor(runs*progress_interval))
                            waitbar(progress, hwaitbar, 'Linear equalization tuning'); figure(hwaitbar); drawnow;
                        end
                    else
                        if ~mod(pxi,floor(runs*progress_interval)), fprintf('%i%% ', round(progress*100) );end
                    end
                    
                    %get the cursor for this iteration
                    txffe_cur=txffe_cursor_vector(TK);
                    
                    % Skip combinations with small values of c(0), not guaranteed to be supported by all transmitters.
                    if txffe_cur<param.tx_ffe_c0_min
                        continue;
                    end
                    old_loops=old_loops+1;
                    
                    %get the index used for each tap on this iteration
                    %this is needed for the LOCAL SEARCH block
                    tx_index_vector=FULL_tx_index_vector(TK,:);
                    
                    %Original LOCAL SEARCH Block:
                    %Keeping this one as commented code because it is a bit more readable than the Modular Block below
                    %But unlike the Modular Block, this one does not work if additional TXFFE taps are added
%                     % speedup "local search" heuristic - Adee Ran 03-17-2020
%                     % skip configurations more than
%                     % 2 steps away from current "best" point on any grid direction
%                     %  Matt Brown 11/19/2021 for cp2 and cp3
%                     if param.LOCAL_SEARCH>0 && ~isinf(best_FOM) && ...
%                                ((k_cp2>1 && length(cp3_values)>1 && abs(k_cp3-find(cp3_values==best_txffe(cur+3)))>param.LOCAL_SEARCH) ...
%                             || (k_cp1>1 && length(cp2_values)>1 && abs(k_cp2-find(cp2_values==best_txffe(cur+2)))>param.LOCAL_SEARCH) ...
%                             || (k_cm1>1 && length(cp1_values)>1 && abs(k_cp1-find(cp1_values==best_txffe(cur+1)))>param.LOCAL_SEARCH) ...
%                             || (k_cm2>1 && length(cm1_values)>1 && abs(k_cm1-find(cm1_values==best_txffe(cur-1)))>param.LOCAL_SEARCH) ...
%                             || (k_cm3>1 && length(cm2_values)>1 && abs(k_cm2-find(cm2_values==best_txffe(cur-2)))>param.LOCAL_SEARCH) ...
%                             || (k_cm4>1 && length(cm3_values)>1 && abs(k_cm3-find(cm3_values==best_txffe(cur-3)))>param.LOCAL_SEARCH) ...
%                             || (g_LP_index>1 && length(cm4_values)>1 && abs(k_cm4-find(cm4_values==best_txffe(cur-4)))>param.LOCAL_SEARCH) ...
%                             || (ctle_index>1 && abs(g_LP_index-best_G_high_pass)>param.LOCAL_SEARCH))
% 
%                         continue;
%                     end

                    %Modular LOCAL_SEARCH block:
                    % speedup "local search" heuristic - Adee Ran 03-17-2020
                    % skip configurations more than 2 steps away from current "best" point on any grid direction
                    skip_it=0;
                    if param.LOCAL_SEARCH>0 && ~isinf(best_FOM)
                        %instead of looping across all taps, only loop across
                        %those with length>1 (txffe_sweep_indices).
                        %It saves time since this block is encountered so often
                        for kj=1:num_txffe_sweep_indices
                            kv=txffe_sweep_indices(kj);
                            if kv==1
                                previous_loop_val=g_LP_index;
                            else
                                previous_loop_val=tx_index_vector(kv-1);
                            end
                            if previous_loop_val>1
                                best_index_this_tap=best_txffe_index(kv);
                                if abs(tx_index_vector(kv)-best_index_this_tap)>param.LOCAL_SEARCH
                                    skip_it=1;
                                    break;
                                end
                            end
                        end
                        
                        if ~skip_it && ctle_index>1 && abs(g_LP_index-best_G_high_pass)>param.LOCAL_SEARCH
                            skip_it=1;
                        end
                    end
                    if skip_it
                        continue;
                    end
                    %End Modular LOCAL SEARCH block
                    
                    new_loops=new_loops+1;
                    
                    %fetch txffe for this iteration
                    txffe=txffe_matrix(TK,:);
                    
                    %The phase shift exponentials used in get_xtlk_noise are independent of
                    %everything except number of taps and cursor position
                    %So it can be calculated 1 time here to avoid thousands of re-calcs
                    if ~calc_exp_phase
                        calc_exp_phase=1;
                        for k=1:length(txffe)
                            phase_memory(:,k)=exp(-1j*2*pi*(k-cur).*f/param.fb);
                        end
                        if OP.RxFFE
                            for k=-1*param.RxFFE_cmx:param.RxFFE_cpx
                                phase_memoryRXFFE(:,k+param.RxFFE_cmx+1)=exp(-1j*2*pi*(k+1).*f/param.fb);
                            end
                            phase_memory=[phase_memory phase_memoryRXFFE];
                        end
                    end
                    
                    %% Unequalized Pulse Reponse & circshift for FFE
                    %Perform circshift for FFE only when CTLE is updated.  The FFE_Fast function takes
                    %in the pre-shifted matrix
                    if ctle_response_updated
                        ctle_response_updated=0;
                        num_pre=cur-1;
                        %Another speed up:  the unequalized pulse is also only unique for each CTLE update
                        %Calculating here reduces number of convolutions by thousands
                        if OP.TDMODE
                            pulse_ctle=chdata(1).(ctle_field)(:);
                        else
                            %uneq_pulse=filter(ones(param.samples_per_ui, 1), 1, chdata(1).(ctle_field)(:));
                            %"conv2" is faster than filter. Just need to chop off extra points at the end
                            pulse_ctle=conv2(chdata(1).(ctle_field)(:),ones(param.samples_per_ui, 1));
                            pulse_ctle=pulse_ctle(1:end-param.samples_per_ui+1);
                        end
                        for k=1:length(txffe)
                            pulse_ctle_circshift(:,k)=circshift(pulse_ctle,[(k-1-num_pre)*param.samples_per_ui 0]);
                        end
                    end
                    
                    %% Apply TXFFE to pre-shifted pulse response
                    %[sbr] = FFE( txffe, cur-1, param.samples_per_ui, uneq_pulse);
                    sbr=FFE_Fast(txffe,pulse_ctle_circshift);
                    sbr_from_txffe=sbr;
                    sbr1=sbr;
                    
                    
                    %% Find Sample Location
                    % If RXFFE is included, the sample location will be found again below
                    [cursor_i,no_zero_crossing,sbr_peak_i,zxi]=cursor_sample_index(sbr,param,OP,start_max_idx:end_max_idx);
                    if param.ts_anchor==0
                        %keep MM
                    elseif param.ts_anchor==1
                        %peak sample
                        cursor_i=sbr_peak_i;
                        no_zero_crossing=0;
                    elseif param.ts_anchor==2
                        %max DV
                        possible_cursor=sbr(sbr_peak_i-param.samples_per_ui:sbr_peak_i+param.samples_per_ui);
                        possible_precursor=sbr(sbr_peak_i-2*param.samples_per_ui:sbr_peak_i);
                        [max_diff,d_idx]=max(possible_cursor-possible_precursor);
                        cursor_i=sbr_peak_i-param.samples_per_ui+d_idx-1;
                        no_zero_crossing=0;
                    else
                        error('ts_anchor parameter must be 0, 1, or 2');
                    end
                    if no_zero_crossing
                        continue;
                    end
                    raw_cursor_i=cursor_i;
                    
                    %%%%%%%%%%
                    %%%%%%%%%%
                    %%%%%%%%%%
                    %NEW ITICK LOOP (not indenting everything yet)
                    [~,si]=sort(abs(full_sample_range));
                    best_positive_itick_FOM=-inf;
                    best_negative_itick_FOM=-inf;
                    best_positive_itick_in_loop=[];
                    best_negative_itick_in_loop=[];
                    best_itick_FOM=-inf;
                    best_itick_in_cluster=[];
                    best_cluster=[];
                    
                    %box_search:  take the middle of 5 point windows, then use the best of those to search the rest of that window
                    %middle_search:  start from itick=0 and work outwards in negative and positive direction. stop searching when FOM starts to go down (using LOCAL SEARCH) 
                    box_search=0;
                    middle_search=1;
                    
                    if box_search
                        box_size=5;
                        box_mid=floor(box_size/2);
                        cluster=full_sample_range(1)+box_mid:box_size:full_sample_range(end);
                        CL=length(cluster);
                        loop_range=1:CL+box_mid*2;
                    elseif middle_search
                        loop_range=si;
                    else
                        loop_range=1:length(full_sample_range);
                    end
                    
                    for itickn=loop_range
                        if box_search
                            if itickn<=CL
                                itick=cluster(itickn);
                            else
                                if itickn==CL+1
                                    best_cluster=setdiff([best_itick_in_cluster-box_mid:best_itick_in_cluster+box_mid],best_itick_in_cluster);
                                end
                                if isempty(best_cluster)
                                    continue;
                                end
                                itick=best_cluster(itickn-CL);
                            end
                        else
                            itick=full_sample_range(itickn); 
                        end
                    
                    itick_cases=itick_cases+1;
                    
                    sbr=sbr_from_txffe;
                    cursor_i = raw_cursor_i+itick;
                    
                    %Local Search for +/- itick sweep
                    if middle_search && param.LOCAL_SEARCH>0
                        if itick>=0 && ~isinf(best_positive_itick_FOM) && abs(best_positive_itick_in_loop-itick)>=param.LOCAL_SEARCH
                            itick_skips=itick_skips+1;
                            continue;
                        end
                        if itick<=0 && ~isinf(best_negative_itick_FOM) && abs(best_negative_itick_in_loop-itick)>=param.LOCAL_SEARCH
                            itick_skips=itick_skips+1;
                            continue;
                        end
                    end
                    
                    triple_transit_time = round(sbr_peak_i*2/param.samples_per_ui)+20;
                    if min_number_of_UI_in_response < triple_transit_time
                        min_number_of_UI_in_response = triple_transit_time;
                    end
                    
                    cursor = sbr(cursor_i);
                    
                    %% RXFFE
                    if OP.RxFFE
                        % [ sbr, C]=force(sbr,param,OP,cursor_i);
                        %[ sbr, C]=force(sbr,param,OP,zxi+param.samples_per_ui);
                        %if isrow(sbr), sbr=sbr';end
                        
                        %AJG:  do not return sbr here (run time improvement)
                        %UPDATE:  use cursor_i in RXFFE instead of zero crossing + 1 UI
                        %[ ~, C]=force(sbr,param,OP,zxi+param.samples_per_ui,[],0);
                        % [ ~, C, floating_tap_locations]=force(sbr,param,OP,cursor_i,[],0);
                        % sbr at this point include the current setting
                        %                under consideration of txffe h21 ctf and fr
                        switch upper(OP.FFE_OPT_METHOD)
                            case 'MMSE'
                                OP.WO_TXFFE=0;
                                PSD_results=get_PSDs(PSD_results,sbr,cursor_i,txffe,gdc_values(ctle_index),g_DC_low,param,chdata,OP);
                                S_n=PSD_results.S_rn+PSD_results.S_tn+PSD_results.S_xn+PSD_results.S_jn;
                                if 0 % for debug
                                    plot(PSD_results.fvec(1:param.num_ui_RXFF_noise)/param.fb,10*log10(PSD_results.S_rn*1000/100) ,'disp','Srn')
                                    hold on
                                    plot(PSD_results.fvec(1:param.num_ui_RXFF_noise)/param.fb,10*log10(PSD_results.S_xn*1000/100) ,'disp','Sxn')
                                    plot(PSD_results.fvec(1:param.num_ui_RXFF_noise)/param.fb,10*log10(PSD_results.S_tn*1000/100) ,'disp','Stn')
                                    plot(PSD_results.fvec(1:param.num_ui_RXFF_noise)/param.fb,10*log10(PSD_results.S_jn*1000/100) ,'disp','Sjn')
                                    plot(PSD_results.fvec(1:param.num_ui_RXFF_noise)/param.fb,10*log10(PSD_results.S_n*1000/100) ,'disp','Sn')
                                     xlim([0 0.5])
                                    % ylim([-190 -160])
                                    set(gcf,'defaulttextinterpreter','none')
                                    xlabel('Normalized Frequency')
                                    ylabel('PSD dBm/Hz')
                                    hold on
                                    grid on
                                    legend show
                                    title('PSD')
                                end
                                MMSE_results = MMSE(PSD_results,sbr,cursor_i, param, OP ) ;
                                % floating_tap_locations=MMSE_results.MLSE_results;
                                C=MMSE_results.C;
                                FOM=MMSE_results.FOM;
                                floating_tap_locations=MMSE_results.floating_tap_locations;
                            otherwise
                                [ ~, C, floating_tap_locations]=force(sbr,param,OP,cursor_i,[],0, chdata, txffe, Noise_XC);
                        end
                        %Now there is a stand alone function for determining if RXFFE taps are illegal
                        %This is because the "force" function will also do a legality check when "backoff" is enabled
                        if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
                            if RXFFE_Illegal(C,param)
                                continue;
                            end
                        end
                        %AJG:  speed up:  calculate sbr after checks for illegal taps (many tap combinations are illegal and "FFE" is time consuming)
                        sbr=FFE(C,param.RxFFE_cmx,param.samples_per_ui,sbr);
                        if isrow(sbr), sbr=sbr';end
                        
                        %% second guess at cursor location (t_s)  - based on approximate zero crossing
                        %This entire block is now inside the "if OP.RxFFE" to avoid unnecessary re-calculation
                        %UPDATE:  NOT RESAMPLING AFTER RXFFE
%                         [cursor_i,no_zero_crossing,sbr_peak_i]=cursor_sample_index(sbr,param,OP,start_max_idx:end_max_idx);
%                         if no_zero_crossing
%                             continue;
%                         end
                        
                        cursor = sbr(cursor_i);
                    end
                    A_p=sbr(sbr_peak_i);
                    %% 93A.1.6 step c defines A_s %%
                    A_s = param.R_LM*cursor/(param.levels-1);
                    if isempty(delta_sbr)
                        delta_sbr = sbr;
                    end
                    sbr=sbr(:);
                    %% Equation 93A-27 "otherwise" case %% param.N_bmax is param.ndfe if groups are not used
                    
                    if(param.Floating_DFE), param.ndfe=param.N_bmax; end
                    far_cursors = sbr(cursor_i-T_O+param.samples_per_ui*(param.ndfe+1):param.samples_per_ui:end);
                    t=((cursor_i+param.samples_per_ui*(param.ndfe+1):param.samples_per_ui:length(sbr))-(cursor_i+param.samples_per_ui*(param.ndfe+1)))*...
                        param.ui/param.samples_per_ui;
                    precursors = sbr(cursor_i-param.samples_per_ui:-param.samples_per_ui:1);
                    precursors = precursors(end:-1:1);
                    
                    %                                 % Error message if the sbr is not long enough for the specified range of Nb
                    %                                 if length(sbr) < cursor_i+param.samples_per_ui*(param.ndfe+1)
                    %                                     close(hwaitbar);
                    %                                     error('Pulse Response contains %d samples after the cursor. Specified Nb requires %d samples after the cursor.' ...
                    %                                         , length(sbr)-cursor_i, param.samples_per_ui*(param.ndfe+1));
                    %                                 end
                    
                    
                    
                    %% skip this case if FOM has no chance of beating old FOM
                    %this is also done below but with excess_dfe_cursors included.
                    %excess_dfe_cursors requires the floating DFE computation which is
                    %time consuming, so checking here can have significant run time improvements
                    sigma_ISI_ignoreDFE = param.sigma_X*norm([precursors;  far_cursors]);
                    if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
                        if (20*log10(A_s/sigma_ISI_ignoreDFE) < best_FOM)
                            continue
                        end
                    end
                    
                    %% Equation 93A-27, when 1<=n<=N_b
                    %required length = cursor + all DFE UI + 1 additional UI
                    sbr_required_length=cursor_i+param.samples_per_ui*(param.ndfe+1);
                    if length(sbr)<sbr_required_length
                        sbr(end+1:sbr_required_length)=0;
                    end
                    dfecursors=sbr(cursor_i+param.samples_per_ui*(1):param.samples_per_ui:cursor_i+param.samples_per_ui*(param.ndfe));
                    if param.dfe_delta ~= 0
                        dfecursors_q=floor(abs(dfecursors/sbr(cursor_i))./param.dfe_delta).*param.dfe_delta.*sign(dfecursors)*sbr(cursor_i);
                        
                    else
                        dfecursors_q=dfecursors;
                    end
                    if param.Floating_DFE
                        %% floating taps
                        postcurors= sbr(cursor_i+param.samples_per_ui:param.samples_per_ui:end);
                        
                        [floating_tap_locations, floating_tap_coef, hisi, bmax]= floatingDFE( postcurors ,param.ndfe_passed,param.N_bf,param.N_bg,param.N_bmax, param.bmaxg, sbr(cursor_i), param.dfe_delta  );
                        
                        newbmax= [ param.bmax bmax(param.ndfe_passed+1:param.N_bmax)].';
                        param.use_bmax=newbmax;
                        %AJG021820
                        param.use_bmin=[param.bmin bmax(param.ndfe_passed+1:param.N_bmax)*-1].';
                    else
                        param.use_bmax=param.bmax;
                        %AJG021820
                        param.use_bmin=param.bmin;
                    end
                    
                    %AJG021820
                    actual_dfecursors=dfe_clipper(dfecursors_q,sbr(cursor_i)*param.use_bmax(:),sbr(cursor_i)*param.use_bmin(:));
                    if do_C2M
                        dfecursors_windowed=sbr(cursor_i-T_O+param.samples_per_ui*(1):param.samples_per_ui:cursor_i+param.samples_per_ui*(param.ndfe)-T_O);
                        % readjust SBR
                        if 0
                            %PR_DFE_center not currently used, so this is in "if 0" statement
                            PR_DFE_center=sbr;
                            for n=1:param.ndfe
                                %                                                 for ix=-param.samples_per_ui/2: param.samples_per_ui/2
                                %                                                     i_sample=ix+n*param.samples_per_ui+cursor_i;
                                %                                                     dper=sbr(i_sample)- actual_dfecursors(n);
                                %                                                     PR_DFE_center(i_sample)=dper;
                                %                                                 end
                                i_sample=(-param.samples_per_ui/2: param.samples_per_ui/2)+n*param.samples_per_ui+cursor_i;
                                PR_DFE_center(i_sample)=sbr(i_sample)-actual_dfecursors(n);
                            end
                        end
                        excess_dfe_cursors=dfecursors_windowed-actual_dfecursors;
                    else
                        excess_dfe_cursors=dfecursors-actual_dfecursors;
                    end
                    dfetaps=actual_dfecursors/sbr(cursor_i);
                    
                    if length(dfetaps) >= param.N_tail_start && param.N_tail_start ~=0
                        tail_RSS=norm(dfetaps(param.N_tail_start:end));
                        if  tail_RSS ~= 0
                            if tail_RSS >= param.B_float_RSS_MAX
                                param.use_bmax(param.N_tail_start:end)= ...
                                    min(tail_RSS, param.B_float_RSS_MAX) * sign(dfetaps(param.N_tail_start:end)).*dfetaps(param.N_tail_start:end) /tail_RSS;
                                %AJG021820
                                param.use_bmin(param.N_tail_start:end)= ...
                                    min(tail_RSS, param.B_float_RSS_MAX) * -1 * sign(dfetaps(param.N_tail_start:end)).*dfetaps(param.N_tail_start:end) /tail_RSS;
                            end
                        end
                        
                        %AJG021820
                        actual_dfecursors=dfe_clipper(dfecursors_q,sbr(cursor_i)*param.use_bmax(:),sbr(cursor_i)*param.use_bmin(:));
                        if do_C2M
                            excess_dfe_cursors=dfecursors_windowed-actual_dfecursors;
                        else
                            excess_dfe_cursors=dfecursors-actual_dfecursors;
                        end
                        dfetaps=actual_dfecursors/sbr(cursor_i);
                        
                    else
                        tail_RSS=0;
                    end
                    %% Eq. 93A-28 %%
                    sampling_offset = mod(cursor_i, param.samples_per_ui);
                    %ensure we can take early sample
                    if sampling_offset<=1
                        sampling_offset=sampling_offset+param.samples_per_ui;
                    end
                    if (OP.LIMIT_JITTER_CONTRIB_TO_DFE_SPAN)
                        cursors_early_sample = sbr(cursor_i-1+param.samples_per_ui*(-1:param.ndfe));
                        cursors_late_sample = sbr(cursor_i+1+param.samples_per_ui*(-1:param.ndfe));
                    else
                        cursors_early_sample = sbr(sampling_offset-1:param.samples_per_ui:end);
                        cursors_late_sample = sbr(sampling_offset+1:param.samples_per_ui:end);
                    end
                    % ensure lengths are equal
                    cursors_early_sample = cursors_early_sample(1:length(cursors_late_sample));
                    h_J = (cursors_late_sample-cursors_early_sample)/2*param.samples_per_ui;
                    if ~OP.SNR_TXwC0
                        %% Equation 93A-30 %%
                        % since A_s = param.R_LM*cursor/(param.levels-1), cursor=(param.levels-1)*A_s/param.R_LM
                        sigma_TX = (param.levels-1)*A_s/param.R_LM*10^(-param.SNR_TX/20);
                    else
                        sigma_TX = (param.levels-1)*A_s/txffe(cur)/param.R_LM*10^(-param.SNR_TX/20);% SNER_TX mod from Adee
                    end
                    %% Equation 93A-31 %%
                    sigma_ISI = param.sigma_X*norm([precursors; excess_dfe_cursors; far_cursors]);
                    ISI_N=param.sigma_X*norm( far_cursors);
                    %% break if FOM has no chance of beating old e
                    OP.exe_mode=1;
                    if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
                        switch OP.EXE_MODE
                            case 0
                            case 1
                                if (20*log10(A_s/sigma_ISI) < best_FOM)
                                    continue
                                end
                            case 2
                                if (20*log10(A_s/sigma_ISI) < best_FOM)
                                    break
                                end
                        end
                    end
                    %% Equation 93A-32 %%
                    sigma_J = norm([param.A_DD param.sigma_RJ])*param.sigma_X*norm(h_J);
                    
                    %% Equations 93A-33 and 93A-34 for FEXT (depends on TXFFE setting) %%
                    if OP.RX_CALIBRATION
                        sigma_XT=0;
                    else
                        if ~OP.RxFFE
                            % sigma_FEXT =  get_xtlk_noise( 0, 'FEXT', param ,chdata);
                            % sigma_XT = norm([sigma_NEXT sigma_FEXT]);
                            [sigma_XT,~,~] =  get_xtlk_noise( txffe, 'FEXT', param ,chdata,phase_memory); %with three outputs, the sigma_XT includes both FEXT and NEXT zhilei huang 01/11/2019
                            %% Equation 93A-36 denominator (actually its sqrt)
                        else % John Ewen: 13/12/20018
                            [sigma_XT,~,~] =  get_xtlk_noise( txffe, 'FEXT', param ,chdata, phase_memory,C); %with three outputs, the sigma_XT includes both FEXT and NEXT zhilei huang 01/11/2019
                            %                                         sigma_FEXT_ffe=  get_xtlk_noise( 0, 'FEXT', param, chdata, C );
                            %                                         if ~OP.RxFFE
                            %                                             sigma_NEXT_ffe =  get_xtlk_noise(0, 'NEXT', param, chdata);
                            %                                         else
                            %                                             sigma_NEXT_ffe =  get_xtlk_noise(0, 'NEXT', param, chdata , C);
                            %                                         end
                            %                                         sigma_XT = norm([sigma_NEXT_ffe sigma_FEXT_ffe])*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
                        end
                    end
                    if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
                        if OP.RxFFE % modify sigma_N with rx noise from the rx ffe
                            index_f2=find(chdata(1).faxis(:)>param.fb,1,'first');
                            if isempty(index_f2), index_f2=length(chdata(1).faxis);end
                            f=chdata(1).faxis;
                            H_Rx_FFE=zeros(1,length(f));
                            for ii=-param.RxFFE_cmx:param.RxFFE_cpx
                                %H_Rx_FFE=C(ii+param.RxFFE_cmx+1).*exp(-1j*2*pi*(ii+1).*f/param.fb)+H_Rx_FFE;
                                if C(ii+param.RxFFE_cmx+1)==0
                                    %speed up:  skip cases when rxffe=0
                                    continue;
                                end
                                if ii+1==0
                                    %speed up:  ii+1=0, so just scalar addition and avoid exp calc
                                    H_Rx_FFE = H_Rx_FFE + C(ii+param.RxFFE_cmx+1);
                                else
                                    H_Rx_FFE=C(ii+param.RxFFE_cmx+1).*transpose(phase_memory(:,ii+param.RxFFE_cmx+1+length(txffe)))+H_Rx_FFE;
                                end
                            end
                            sigma_N =sqrt(param.eta_0*sum(H_sy(2:end) .* abs(H_r(2:end) .* H_ctf(2:end) .*H_Rx_FFE(2:end)).^2 .* diff(chdata(1).faxis)/1e9)); % changed from /chdata(1).faxis(end) B. Kirkland S. Elnagar 11/6/2021
                        end
                    end
                    %% Equation 93A-36 (note log argument is voltage rather than power ratio)
                    total_noise_rms  = norm([sigma_ISI sigma_J sigma_XT sigma_N sigma_TX sigma_ne]);
                    if do_C2M
                        if param.Noise_Crest_Factor == 0
                            ber_q = sqrt(2)*erfcinv(2*param.specBER);
                        else
                            ber_q=param.Noise_Crest_Factor;
                        end
                        if OP.force_pdf_bin_size
                            delta_y = OP.BinSize;
                        else
                            delta_y = min(A_s/1000, OP.BinSize);
                        end
                        ne_noise_pdf = normal_dist(0, ber_q, delta_y);
                        cci_pdf = normal_dist(0, ber_q, delta_y);
                        chdata(1).eq_pulse_response=sbr;
                        tmp_result.t_s= cursor_i;
                        tmp_result.A_s=A_s;
                        EH_1st= 2*(A_s-erfcinv(param.specBER*2)*2/sqrt(2)*total_noise_rms);
                        if EH_1st <=  param.Min_VEO_Test/1000 -.001
                            %                                         sprintf( '%g.1 As ..  %g.1 EH\n',A_s*1000,EH_1st*1000)
                            continue
                        else
                            %                                          sprintf( ' OK before %g As ..  %g EH\n',A_s*1000,EH_1st*1000)
                        end
                        Struct_Noise.sigma_N=sigma_N;
                        Struct_Noise.sigma_TX=sigma_TX;
                        Struct_Noise.cci_pdf=cci_pdf;
                        Struct_Noise.ber_q=ber_q;
                        Struct_Noise.ne_noise_pdf=ne_noise_pdf;
                        [Left_EW,Right_EW,eye_contour,EH_T_C2M,EH_B_C2M]=COM_eye_width(chdata,delta_y,tmp_result,param,OP,Struct_Noise,1);
                        EH=EH_T_C2M-EH_B_C2M;
                        N_i=(A_s*2-EH)/2;
                        if EH <= param.Min_VEO_Test/1000
                            %                                         sprintf( 'After  As=%.1f ..  EH=%.1f  EH_1st=%.1f  \n',A_s*1000,EH*1000, EH_1st*1000)
                            continue
                        else
                            %                                         sprintf( '<strong> After  As=%.1f ..  EH=%.1f  EH_1st=%.1f </strong> \n',A_s*1000,EH*1000, EH_1st*1000)
                        end
                        if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE) % MMSE defines its own FOM
                            FOM =20*log10(A_s/N_i);
                        end
                    else
                        if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE) % MMSE defines its own FOM
                            FOM = 20*log10(A_s/total_noise_rms);
                        end
                        %                                             if strfind(param.CTLE_type,'CL120e')
                        %                                                 FOM = A_s*C(param.RxFFE_cmx+1)-total_noise_rms_ffe;
                    end
                    if 0 % for loop analysis
                        result.FOM_array(new_loops)=FOM;
                    end
                    
                    if FOM>best_itick_FOM
                        best_itick_FOM=FOM;
                        best_itick_in_cluster=itick;
                    end
                    
                    if itick>=0 && FOM>best_positive_itick_FOM
                        best_positive_itick_FOM=FOM;
                        best_positive_itick_in_loop=itick;
                    end
                    if itick<=0 && FOM>best_negative_itick_FOM
                        best_negative_itick_FOM=FOM;
                        best_negative_itick_in_loop=itick;
                    end
                    
                    itick_index=find(itick==full_sample_range);
                    FOM_TRACKER(Gffe_index,ctle_index,g_LP_index,TK,itick_index)=FOM;
                    
                    if (FOM > best_FOM)
                        best_current_ffegain=param.current_ffegain;
                        best_txffe = txffe;
                        %along with best_txffe, save the indices of the best_txffe
                        %(saves time in LOCAL SEARCH block)
                        best_txffe_index=tx_index_vector;
                        best_sbr = sbr;
                        best_ctle = ctle_index;
                        best_G_high_pass =g_LP_index;
                        best_FOM = FOM;
                        best_cursor_i = cursor_i;
                        best_itick = itick;                       
                        if ~OP.TDMODE
                            [ effective_channel ] = FFE( txffe , cur-1, param.samples_per_ui, chdata(1).ctle_imp_response );
                            best_IR=effective_channel;
                        end
                        best_sigma_N = sigma_N;
                        best_h_J = h_J;
                        best_A_s=A_s;
                        best_A_p=A_p;
                        best_ISI=ISI_N;
                        best_bmax=param.use_bmax;
                        %AJG021820
                        best_bmin=param.use_bmin;
                        best_tail_RSS=tail_RSS;
                        best_dfetaps=dfetaps;
                        if param.Floating_DFE
                            best_floating_tap_locations=floating_tap_locations;
                            best_floating_tap_coef=floating_tap_coef;
                        end
                        if param.Floating_RXFFE
                            best_floating_tap_locations=floating_tap_locations;
                            % best_floating_tap_coef=floating_tap_coef;
                        end                        
                        if OP.RxFFE
                            best_RxFFE=C;
                            best_PSD_results=PSD_results;
                            best_MMSE_results=MMSE_results;
                        end
                    end
                    end
                end
                
            end
        end
    end
    if do_C2M
        if  best_FOM == -inf
            param.Min_VEO_Test=0;
        else
            break
        end
    end
end
if 0
    fprintf('old loops = %d\n',old_loops);
    fprintf('new loops = %d\n',new_loops);
    display(sprintf('\n :loops = %g',pxi))
end

%turn this on to review if FOM changes sign more than once in an itick loop
if 0
    DIR_CHANGE={};
    for m=1:length(Gffe_values)
        for n=1:length(gdc_values)
            for k=1:lf_indx
                FOM_this_mat=squeeze(FOM_TRACKER(m,n,k,:,:));
                %x reveals if FOM on a particular row (locked txffe, moving itick) goes up or down
                %1 = goes up, -1=goes down
                x=sign(diff(FOM_this_mat')');
                %y = change in sign on x.  the location of a "2" is where FOM changes direction
                y=abs(diff(x'))';
                %the goal is the FOM only changes direction once. so count the occurences of the 2
                for j=1:size(FOM_this_mat,1)
                    z{j}=find(y(j,:)==2);
                end
                zL=cellfun('length',z);
                %return any row where FOM changed direction more than once
                DIR_CHANGE{j,k}=find(zL>1);
            end
        end
    end
    multi_direction_change=find(~cellfun('isempty',DIR_CHANGE))
end

if ~exist('best_cursor_i', 'var')% take last setting
    result.eq_failed=true;
    display('equalization failed')
    best_bmax=param.bmax;
    %AJG021820
    best_bmin=param.bmin;
    best_tail_RSS=0;
    best_current_ffegain=0;
    best_txffe = txffe;
    best_sbr = sbr;
    best_ctle = ctle_index;
    if OP.RxFFE
        best_PSD_results=PSD_results;
        best_MMSE_results=MMSE_results;
        best_RxFFE=C;
    end
    best_G_high_pass =g_LP_index;
    best_FOM = FOM;
    %if this block is reached, the last encountered EQ parameters are used
    %if it so happened that there was no zero crossing in the last encountered EQ set, then cursor_i will be empty
    %EQ search has failed, so it is not important to give an exact sample location, so just use the peak of the pulse
    if isempty(cursor_i)
        [~,cursor_i]=max(sbr);
    end
    best_cursor_i = cursor_i;
    best_itick = itick;
    if ~OP.TDMODE
        [ effective_channel ] = FFE( txffe , cur-1, param.samples_per_ui, chdata(1).ctle_imp_response );
        best_IR=effective_channel;
    end
    best_sigma_N = sigma_N;
    best_h_J = h_J;
    best_A_p=max(sbr);
    best_ISI=1;
    best_dfetaps= sbr( cursor_i+param.samples_per_ui*(1):param.samples_per_ui:cursor_i+param.samples_per_ui*(param.ndfe))/sbr(cursor_i) ;
    best_A_s= sbr( cursor_i);
    if param.Floating_DFE
        best_floating_tap_locations=[];
        best_floating_tap_coef=[];
    end
    if do_C2M
        return
    end
    %     return
else
    result.eq_failed=false; % RIM 12/30/2023
end

best_cursor = best_sbr(best_cursor_i);
% report during debug
PRin=filter(ones(param.samples_per_ui, 1),1, chdata(1).uneq_imp_response);
%If sbr was zero padded, then PRin needs to do so as well)
if length(PRin)<length(best_sbr)
    PRin(end+1:length(best_sbr))=0;
end
f=1e8:1e8:100e9;

H_bt=Bessel_Thomson_Filter(param,f,OP.Bessel_Thomson);
H_bw=Butterworth_Filter(param,f,OP.Butterworth);
H_RCos=Raised_Cosine_Filter(param,f,OP.Raised_Cosine);% experiment with RCos
% need to include H_RCos in noise and when computing the system ir for thru
% and crosstalk
H_r=H_bw.*H_bt.*H_RCos;

ctle_gain1 = (10^(gdc_values(best_ctle)/20) + 1i*f/param.CTLE_fz(best_ctle)) ./ ...
    ((1+1i*f/param.CTLE_fp1(best_ctle)).*(1+1i*f/param.CTLE_fp2(best_ctle)));

switch param.CTLE_type
    case 'CL93'
        H_low=1;
    case 'CL120d'
        H_low=(10^(param.g_DC_HP_values(best_G_high_pass)/20) +  1i*f/param.f_HP(best_G_high_pass))./(1 + 1i*f/param.f_HP(best_G_high_pass));
    case 'CL120e'
        H_low=(1 +  1i*f/param.f_HP_Z(best_ctle))./ (1 + 1i*f/param.f_HP_P(best_ctle));
end
ctle_gain=H_low.*ctle_gain1.*H_r;



%lsbr=length(sbr);
%use length of best_sbr in case zero padding was performed
%check "sbr_required_length" variable
lsbr=length(best_sbr);
t=0:param.ui/param.samples_per_ui:(lsbr-1)*param.ui/param.samples_per_ui;

sampled_best_sbr_precursors_t = (best_cursor_i/param.samples_per_ui:-1:1/param.samples_per_ui)*param.ui;
sampled_best_sbr_precursors_t = sampled_best_sbr_precursors_t(end:-1:2); % exclude cursor
sampled_best_sbr_precursors = best_sbr(round(sampled_best_sbr_precursors_t/param.ui*param.samples_per_ui));
sampled_best_sbr_postcursors_t = (best_cursor_i:param.samples_per_ui:lsbr)/param.samples_per_ui*param.ui;
sampled_best_sbr_postcursors_t = sampled_best_sbr_postcursors_t(2:end); % exclude cursor
sampled_best_sbr_postcursors = best_sbr(round(sampled_best_sbr_postcursors_t/param.ui*param.samples_per_ui));
sampled_best_sbr_dfecursors_t = (best_cursor_i/param.samples_per_ui+(1:param.ndfe_passed))*param.ui;
if param.Floating_DFE
    sampled_best_sbr_fdfecursors_t = (best_cursor_i/param.samples_per_ui+(best_floating_tap_locations))*param.ui;
end
% apply max tap value constraint
dfe_cursors =    sampled_best_sbr_postcursors(1:param.ndfe);
dfe_SBRcursors = sampled_best_sbr_postcursors(1:param.ndfe);
if  isrow(best_bmax) == 1, best_bmax=best_bmax.';end

%AJG021820
if  isrow(best_bmin) == 1, best_bmin=best_bmin.';end
DFE_taps_mV=dfe_clipper(dfe_cursors,best_cursor*best_bmax(1:param.ndfe),best_cursor*best_bmin(1:param.ndfe));
if param.Floating_DFE
    FDFE_taps_mV=DFE_taps_mV(best_floating_tap_locations);
end

sampled_best_sbr_postcursors(1:param.ndfe) = dfe_SBRcursors-DFE_taps_mV;
Symbol_Adj = (param.levels-1);% 3A.1.6
if OP.DEBUG ~=0
    if OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
        % display pulse responses in one axis per test case.
        switch upper(OP.TIME_AXIS)
            case 'S' % RIM 11-13-2023 added user selectable xaxis
                xnorm=1;
                xaxis_label='seconds';
                offset=0;
            case 'UI'
                xnorm=param.ui;
                xaxis_label='UI';
                offset=t(best_cursor_i)/xnorm;
            otherwise
                xnorm=1;
                xaxis_label='seconds';
                offset=0;      
        end
        figure_name = sprintf('PKG %d: Equalization effect: %s : ',  OP.pkg_len_select( param.package_testcase_i), param.base);
        fig=findobj('Name', figure_name);
        if isempty(fig), fig=figure('Name', figure_name); end
        figure(fig);set(gcf,'Tag','COM');
        movegui(fig,'north')
        %figure(fig.Number);
        % hax = subplot(length(OP.pkg_len_select), 1, param.package_testcase_i);
        if OP.RxFFE
            ax1=subplot(2,1,1);
        end
        plot(t/xnorm-offset,best_sbr/Symbol_Adj,'disp', '1/2 Symbol 0-3 Equalized PR');
        hold on
        
        PRplt(1:param.samples_per_ui+5)=PRin(1); % line up with the ffe introduced delay
        PRplt(param.samples_per_ui+6:length(t))=PRin(1:length(t)-param.samples_per_ui-5);
        plot((t-param.ui-t(6))/xnorm-offset,PRplt/Symbol_Adj,'r','disp', '1/2 Symbol 0-3 Unequalized (tp5d) PR');
        stem(t(best_cursor_i)/xnorm-offset,best_sbr(best_cursor_i)/Symbol_Adj,'g','disp','Cursor (sample point)');
        title(sprintf('PKG Case %d',  OP.pkg_len_select( param.package_testcase_i)));
        ylabel('volts')
        xlabel(xaxis_label)
        grid on
        legend show
        legend( 'Location', 'best')
        plot((sampled_best_sbr_precursors_t-param.ui/param.samples_per_ui)/xnorm-offset, sampled_best_sbr_precursors'/Symbol_Adj, 'kx', 'disp','Pre cursors');
        plot((sampled_best_sbr_postcursors_t-param.ui/param.samples_per_ui)/xnorm-offset, sampled_best_sbr_postcursors'/Symbol_Adj, 'ko', 'disp','Post cursors');
        if param.ndfe_passed ~=0
            stem((sampled_best_sbr_dfecursors_t-param.ui/param.samples_per_ui)/xnorm-offset,DFE_taps_mV(1:param.ndfe_passed)/Symbol_Adj','m', 'LineWidth',2,'disp','DFE-canceled cursors');
        end
        if param.Floating_DFE
            stem((sampled_best_sbr_fdfecursors_t-param.ui/param.samples_per_ui)/xnorm-offset,FDFE_taps_mV/Symbol_Adj, 'MarkerFaceColor','red','MarkerEdgeColor','m','LineWidth',1,'disp','FDFE-canceled cursors');
        end
        if OP.RxFFE
            ax2=subplot(2,1,2);
            if param.Floating_RXFFE
                stem((t(best_cursor_i+(best_floating_tap_locations)*param.samples_per_ui))/xnorm-offset,best_RxFFE(best_floating_tap_locations)...
                    ,'filled','disp','RxFFE floating FFE taps')
                hold on
            end
            stem((t(best_cursor_i+param.samples_per_ui*(-param.RxFFE_cmx:param.RxFFE_cpx)))/xnorm-offset,best_RxFFE(1:param.RxFFE_cmx+param.RxFFE_cpx+1)...
                ,'filled','disp','RxFFE fixted FFE taps')
            legend show
            zoom xon
            linkaxes([ax1 ax2],'x')
        end


        grid on
        legend show
        legend( 'Location', 'best')
        zoom xon
        % set(hax, 'tag', 'EQE');
        %
        figure(110);set(gcf,'Tag','COM');
        set(gcf, 'Name', 'CTLE selection');
        movegui(gcf, 'southeast');
        semilogx(f,20*log10(abs(ctle_gain)), 'disp', sprintf('Case %d', param.package_testcase_i));
        hold on
        semilogx(f,20*log10(abs(H_r)), 'disp', sprintf('Rx filter Case %d', param.package_testcase_i));
        semilogx(f,20*log10(abs(H_low.*ctle_gain1)), 'disp', sprintf('CTF Case %d', param.package_testcase_i));
        fbaud_tick=find(f >= baud_rate, 1);
        fnq_tick=find(f >= baud_rate/2, 1);
        stem(f(fnq_tick),20*log10(abs(ctle_gain(fnq_tick))),'g', 'handlevisibility', 'off');
        stem(f(fbaud_tick),20*log10(abs(ctle_gain(fbaud_tick))),'g', 'handlevisibility', 'off');
        recolor_plots(gca);
        title('CTF/w Rx Filter Response')
        ylabel('dB')
        xlabel('Hz')
        legend show
    end
    display(['FOM:                ' ,num2str(best_FOM, 2),' dB']);
    display(['TXFFE coefficients: ' ,mat2str(best_txffe) ] );
    display(['SNR ISI:                ' ,num2str(20*log10(best_A_p/best_ISI), 2),' dB']);
    display(['CTLE DC gain:       ' ,num2str(gdc_values(best_ctle)), ' dB']);
    display(['CTF peaking gain:  ' ,num2str(20*log10(max(abs(ctle_gain))), 2), ' dB']);
    display(['Symbol Available signal:   ' ,num2str(best_cursor/Symbol_Adj)]);
end
if OP.DEBUG && OP.DISPLAY_WINDOW && OP.RX_CALIBRATION==0
    eqe_axes = findobj('tag', 'EQE');
    if ~isempty(eqe_axes), linkaxes(eqe_axes, 'xy'); end
end
if OP.DISPLAY_WINDOW
    close(hwaitbar);
else
    fprintf('\n');
end

% % eq_data
result.cur=cur;
result.txffe = best_txffe;
result.ctle = best_ctle;
result.best_G_high_pass=best_G_high_pass;
result.DFE_taps = best_dfetaps; %relative
result.DFE_taps_i = best_cursor_i+(1:param.ndfe)*param.samples_per_ui;
if param.Floating_DFE
    result.floating_tap_locations=best_floating_tap_locations;
    result.floating_tap_coef=best_floating_tap_coef;
end
if param.Floating_RXFFE
    result.floating_tap_locations=best_floating_tap_locations;
end
result.A_s = best_A_s;
result.t_s = best_cursor_i;
result.itick = best_itick;
result.sigma_N = best_sigma_N;
result.h_J = best_h_J;
result.FOM = best_FOM;
if ~OP.TDMODE
    %If sbr was zero padded, then best_IR needs to do so as well)
    if length(best_IR)<length(best_sbr)
        best_IR(end+1:length(best_sbr))=0;
    end
    result.IR = best_IR;
end
result.t=t;
result.sbr=best_sbr;
if OP.RxFFE
    result.RxFFE=best_RxFFE;
    result.PSD_results=best_PSD_results;
    result.MMSE_results=best_MMSE_results;
end



% changed RIM 4/17/2019 use sum(IR) for Vf of originl IR at N_b UI
% updated RIM 12/17/2021 
result.A_p = max(chdata(1).uneq_pulse_response);
its=find(chdata(1).uneq_pulse_response>=max(chdata(1).uneq_pulse_response),1,'first');
PR=chdata(1).uneq_pulse_response;
iend = its+param.N_v*param.samples_per_ui-param.samples_per_ui/2; % from eq: 163A-3
ibeg= its-param.D_p*param.samples_per_ui-param.samples_per_ui/2;% from eq: 163A-3
if iend >= length(PR)
    iend = length (PR);
end
if ibeg < 1
    ibeg = 1;
end
PR=PR(ibeg:iend);
result.A_f = sum(PR/param.samples_per_ui); %% eq 163A-3
SRn=PR;
for ik=1:floor(length(PR)/param.samples_per_ui)
    SPR=circshift(PR,param.samples_per_ui*ik);
    SPR(1:ik*param.samples_per_ui)=0;
    SRn=SRn+ SPR;
end
codedebug=0;
if codedebug
    fig=figure('Name', 'step and pulse response for code debug');
    figure(fig);set(gcf,'Tag','COM');
    UI=(1:length(SRn))/param.samples_per_ui-param.D_p;
    plot(UI,SRn)
    hold on
    plot(UI,PR)
    xlim([-param.D_p param.N_v])
    grid on;hold off;
    result.step=SRn;
end
i20=find(SRn>=0.20*result.A_f,1,'first');
i80=find(SRn>=0.80*result.A_f,1,'first');
result.Tr_measured_from_step=(i80-i20)/(param.fb*param.samples_per_ui);
result.Pmax_by_Vf=result.A_p/result.A_f;
result.ISI =best_ISI;
result.SNR_ISI=20*log10(best_A_p/best_ISI);
result.best_current_ffegain=best_current_ffegain;
result.best_bmax=best_bmax;
%AJG021820
result.best_bmin=best_bmin;
result.tail_RSS=best_tail_RSS;
function param=parameter_size_adjustment(param,OP)

make_length2={'C_pkg_board' 'C_diepad' 'L_comp' 'C_bump' 'tfx' 'C_v' 'C_0' 'C_1' 'pkg_Z_c' 'brd_Z_c' 'R_diepad'};
make_length_WCPORTZ={'a_thru' 'a_fext' 'a_next' 'SNDR'};
make_length_GDC={'CTLE_fp1' 'CTLE_fp2' 'CTLE_fz' 'f_HP_Z' 'f_HP_P'};
make_length_DCHP={'f_HP'};
make_length_ncases={'AC_CM_RMS'};

%ncases used by make_length_ncases fields
[ncases, mele]=size(param.z_p_tx_cases); % need find the number of test cases RIM 01-08-20

%PORTZ_mult used by make_length_WCPORTZ fields
pkg_sel_vec=ones(1,max(OP.pkg_len_select));
if OP.WC_PORTZ
    PORTZ_mult=[1 1];
else
    PORTZ_mult=pkg_sel_vec;
end

%Parameters that have length = 2
for j=1:length(make_length2)
    if numel(param.(make_length2{j}))==1
        param.(make_length2{j}) = param.(make_length2{j})*[1 1];
    end
end

%Parameters that have length = ncases
for j=1:length(make_length_ncases)
    if numel(param.(make_length_ncases{j}))==1
        param.(make_length_ncases{j}) = param.(make_length_ncases{j})*ones(1,ncases);
    end
end

%Parameters that have length = length(ctle_gdc_values)
for j=1:length(make_length_GDC)
    if numel(param.(make_length_GDC{j}))==1
        param.(make_length_GDC{j}) = param.(make_length_GDC{j})*ones(size(param.ctle_gdc_values));
    end
end

%Parameters that have length = length(g_DC_HP_values)
for j=1:length(make_length_DCHP)
    if numel(param.(make_length_DCHP{j}))==1
        param.(make_length_DCHP{j}) = param.(make_length_DCHP{j})*ones(size(param.g_DC_HP_values));
    end
end

%Parameters that have length associated with PORTZ_mult
for j=1:length(make_length_WCPORTZ)
    if numel(param.(make_length_WCPORTZ{j}))==1
        param.(make_length_WCPORTZ{j}) = param.(make_length_WCPORTZ{j})*PORTZ_mult;
    end
end
function sgm = pdf2sgm(pdf)
avg = sum(pdf.x .* pdf.y);
sgm = sqrt(sum((pdf.x - avg).^2 .* pdf.y));
% end yasuo patch


%% adding tx packgage
function cdf=pdf_to_cdf(pdf)

%Transform PDF to CDF
%The CDF is natively calculated from negative-to-positive voltage.
%This only gives BER calculation for bottom eye.  Need to also
%calculate a CDF of reversed PDF to get top eye.  The final CDF is the
%min of top and bottom CDF values.
%If only interested in one side, a simple cumsum on y is all that is needed.

cdf.yB=cumsum(pdf.y);
cdf.yT=fliplr(cumsum(fliplr(pdf.y)));
cdf.y=min([cdf.yB(:) cdf.yT(:)],[],2);
cdf.x=pdf.x;
function plot_bathtub_curves(hax, max_signal, sci_pdf, cci_pdf, isi_and_xtalk_pdf, noise_pdf,jitt_pdf, combined_interference_and_noise_pdf, bin_size)
cursors = d_cpdf(bin_size,max_signal*[-1 1], [1 1]/2);
signal_and_isi_pdf = conv_fct(cursors, sci_pdf);
signal_and_xtalk_pdf = conv_fct(cursors, cci_pdf);
signal_and_channel_noise_pdf = conv_fct(cursors, isi_and_xtalk_pdf);
signal_and_system_noise_pdf = conv_fct(cursors, noise_pdf);
signal_and_system_jitt_pdf = conv_fct(cursors, jitt_pdf);
signal_and_total_noise_pdf = conv_fct(cursors, combined_interference_and_noise_pdf);
%% Added by Bill Kirkland, June 14, 2017
cursors_l = cursors; cursors_l.y(cursors_l.x>0) = 0;
cursors_r = cursors; cursors_r.y(cursors_r.x<0) = 0;
signal_and_total_noise_pdf_l = conv_fct(cursors_l, combined_interference_and_noise_pdf);
signal_and_total_noise_pdf_r = conv_fct(cursors_r, combined_interference_and_noise_pdf);

semilogy(signal_and_isi_pdf.x, abs(cumsum(signal_and_isi_pdf.y)-0.5) ,'r','Disp','ISI', 'parent', hax)
hold on
semilogy(signal_and_xtalk_pdf.x, abs(cumsum(signal_and_xtalk_pdf.y)-0.5) ,'b','Disp','Xtalk', 'parent', hax)
semilogy(signal_and_channel_noise_pdf.x, abs(cumsum(signal_and_channel_noise_pdf.y)-0.5) ,'c','Disp','ISI+Xtalk', 'parent', hax)
semilogy(signal_and_system_noise_pdf.x, abs(cumsum(signal_and_system_noise_pdf.y)-0.5) ,'m','Disp','Jitter, SNR_TX,RL_M, eta_0 noise', 'parent', hax)
semilogy(signal_and_system_jitt_pdf.x, abs(cumsum(signal_and_system_jitt_pdf.y)-0.5) ,'g','Disp','Jitter noise', 'parent', hax)

%% Added by Bill Kirkland, June 14, 2017
% modification allows bathtub curves to cross over and hence one can
% directly read the noise component.
%semilogy(signal_and_total_noise_pdf.x, abs(cumsum(signal_and_total_noise_pdf.y)-0.5) ,'k','Disp','total noise PDF', 'parent', hax)
vbt_l = abs(0.5-cumsum(signal_and_total_noise_pdf_l.y));
vbt_r = fliplr(0.5-(cumsum(fliplr(signal_and_total_noise_pdf_r.y))));
semilogy(signal_and_total_noise_pdf_l.x, vbt_l ,'k','Disp','total noise PDF left', 'parent', hax)
semilogy(signal_and_total_noise_pdf_r.x, vbt_r ,'k','Disp','total noise PDF right', 'parent', hax)

hc=semilogy(max_signal*[-1 -1 1 1], [0.5 1e-20 1e-20 0.5], '--ok');
set(get(get(hc,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

ylabel(hax, 'Probability')
xlabel(hax, 'volts')
legend(hax, 'show')
% testing code
if 0
    figure_name =  'COM curves';
    fig=findobj('Name', figure_name);
    if isempty(fig), fig=figure('Name', figure_name); end
    figure(fig);set(gcf,'Tag','COM');
    grid on
    semilogy(db(max_signal./(signal_and_total_noise_pdf_r.x)), 2*vbt_r ,'Disp','SNR CDF')
    hold on
    semilogy(db(max_signal./(max_signal-signal_and_total_noise_pdf_r.x)), 2*vbt_r ,'Disp','COM CDF')
    ylim([ 1e-6 0.25])
    xlim([0 30])
    grid on
end
function plot_pie_com(hax, max_signal, sci_pdf, cci_pdf, isi_and_xtalk_pdf, noise_pdf, combined_interference_and_noise_pdf, bin_size,param)
BER=param.specBER;
delta_dB=param.delta_IL;

iex.combined_interference_and_noise_pdf=find(abs(cumsum(combined_interference_and_noise_pdf.y))>= BER, 1, 'first');
iex.noise_pdf=find(abs(cumsum(noise_pdf.y))>= BER, 1, 'first');
iex.cci_pdf=find(abs(cumsum(cci_pdf.y))>= BER, 1, 'first');
iex.sci_pdf=find(abs(cumsum(sci_pdf.y))>= BER, 1, 'first');


maxn(1)=abs(sci_pdf.x(iex.sci_pdf));
maxn(2)=abs(noise_pdf.x(iex.noise_pdf));
maxn(3)=abs(cci_pdf.x(iex.cci_pdf));
maxn_tot=abs(combined_interference_and_noise_pdf.x(iex.combined_interference_and_noise_pdf));

COM=20*log10(max_signal/maxn_tot);
COM_per_noise(1:3)=COM*(maxn(1:3).^2/sum(maxn.^2));
COM_ISI=sprintf(  '%.2g%%',100*COM_per_noise(1)/sum(COM_per_noise));
COM_SYS=sprintf(  '%.2g%%',100*COM_per_noise(2)/sum(COM_per_noise));
COM_xtalk=sprintf('%.2g%%',100*COM_per_noise(3)/sum(COM_per_noise));

pfctr=exp(-0.09054*delta_dB);% less loss
mfctr=exp(0.09054*delta_dB); % more loss

%less loss
plus_maxn(1)=abs(sci_pdf.x(iex.sci_pdf))*pfctr;
plus_maxn(2)=abs(noise_pdf.x(iex.noise_pdf));
plus_maxn(3)=abs(cci_pdf.x(iex.cci_pdf))*pfctr;
% plus_maxn_tot=(maxn(1)*pfctr+maxn(2)+maxn(3)*pfctr)*maxn_tot/sum(plus_maxn);
plus_maxn_tot=norm(plus_maxn);

minus_maxn(1)=abs(sci_pdf.x(iex.sci_pdf))*mfctr;
minus_maxn(2)=abs(noise_pdf.x(iex.noise_pdf));
minus_maxn(3)=abs(cci_pdf.x(iex.cci_pdf))*mfctr;
% minus_maxn_tot=(maxn(1)*mfctr+maxn(2)+maxn(3)*mfctr)*maxn_tot/sum(minus_maxn);
minus_maxn_tot=norm(minus_maxn);

% more loss
COMp=20*log10(max_signal*pfctr/maxn_tot);
COMp_per_noise(1:3)=COMp*(plus_maxn(1:3).^2/plus_maxn_tot^2);
COMp_ISI=sprintf(  '%.2gdB',COMp_per_noise(1));
COMp_SYS=sprintf(  '%.2gdB',COMp_per_noise(2));
COMp_xtalk=sprintf('%.2gdB',COMp_per_noise(3));
% less loss
COMm=20*log10(max_signal*mfctr/maxn_tot);
COMm_per_noise(1:3)=COMm*(minus_maxn(1:3).^2/minus_maxn_tot^2);
COMm_ISI=sprintf(  '%.2gdB',COMm_per_noise(1));
COMm_SYS=sprintf(  '%.2gdB',COMm_per_noise(2));
COMm_xtalk=sprintf('%.2gdB',COMm_per_noise(3));

xax=[[delta_dB delta_dB delta_dB];[ 0 0 0 ]; [ -delta_dB -delta_dB -delta_dB]];


if(COM<0)
    return
end

labels= { ['ISI ' COM_ISI] ['System noise/jitter ' COM_SYS] [' Crosstalk ' COM_xtalk]};

% pie(COM_per_noise,labels)
% legend(labels,{'ISI COM','System noise/jitter COM','Crosstalk COM'});
% legend('show','Location','bestoutside')
nullbar= [ 0 0 0 ];
bar(xax,[nullbar ;COM_per_noise; nullbar],'stacked')
% bar(xax,[COMp_per_noise;COM_per_noise;COMm_per_noise],'stacked')
hold on
bar(xax,[[COMp 0 0 ];nullbar;nullbar],'stacked','w','barwidth',.5)
bar(xax,[nullbar ;nullbar; [0 0 COMm]],'stacked','w','barwidth',.5)

plot([-delta_dB,0, delta_dB], [param.pass_threshold ,param.pass_threshold ,param.pass_threshold ],'r')
% ax=gca;
% ax.XTickLabels = {'decreasing loss','-','increasing loss'};
set(gca,'XTickLabel', {'decreasing loss','-','increasing loss'})
grid on
legend(labels,'Location','north')
xlabel(['if loss is changed by ' sprintf('%.2gdB',delta_dB) ])
ylabel('COM (dB)')
hold off




function [chdata, param] = process_sxp(param, OP, chdata, SDDch)
num_files=length(chdata);
if ~OP.DISPLAY_WINDOW, fprintf('reading file '); end
for i=1:num_files
    if param.package_testcase_i==1  && i==1
        if OP.TDR && i==1
            S.Frequencies=chdata(i).faxis;
            S.Impedance=100;
            if ~OP.SHOW_BRD
                Sfield='_orig';
            else
                Sfield='_raw';
            end
            S.Parameters(1,1,:)=chdata(i).(['sdd11' Sfield]) ;
            if ~param.FLAG.S2P
                S.Parameters(1,2,:)=chdata(i).(['sdd12' Sfield]) ;
                S.Parameters(2,1,:)=chdata(i).(['sdd21' Sfield]) ;
                S.Parameters(2,2,:)=chdata(i).(['sdd22' Sfield]) ;
                S.NumPorts=2; % rim 2/26/2019 correct from S.NumPorts=4;
            else
                S.NumPorts=1;
            end
            if OP.TDR_W_TXPKG
                if OP.ERL == 2
                    error('Cannot add pacakge to s2p files. ERL==2 not supportted if TDR_W_TXPKG = 1')
                end
                R_diepad = param.R_diepad;
                % RX package length is assumed to be the same for all
                % channel types.  swap Cp and Cd for Tx. TDR_W_TXPKG only
                % for  Rx pkg
                [ s11in, s12in, s21in, s22in]=make_full_pkg('TX',S.Frequencies,param,'THRU');
                [ S.Parameters(1,1,:), S.Parameters(1,2,:), S.Parameters(2,1,:), S.Parameters(2,2,:)]  = ...
                    combines4p(  s11in, s12in, s21in, s22in, ...
                    S.Parameters(1,1,:), S.Parameters(1,2,:), S.Parameters(2,1,:), S.Parameters(2,2,:) );
                %                     S=sparameters(S.Parameters,S.Frequencies,100);
                S=SL(S,S.Frequencies,R_diepad(1)*2);
                chdata(i).TX_RL=S.Parameters(2,2,:);
                S.Parameters(1,1,:)=chdata(i).sdd11_orig; % when looking at Tx don't include package
            end
            
            % need to combine S wiht is page and channel
            if param.FLAG.S2P
                port_sel=1;
            else
                port_sel=[1 2];
                if OP.AUTO_TFX
                    [ fir4del, tu] =get_RAW_FIR(squeeze(chdata(i).sdd12_orig),S.Frequencies,OP,param);
                    pix=find(fir4del==max(fir4del),1);
                     param.tfx(2)=2*tu(pix);
                end
            end
            OP.impulse_response_truncation_threshold=1e-5; %Only for TDR not returned out of "process_sxp" function
            for ipsl=1:length(port_sel) % do for both port if s4p
                for izt=1:length(param.Z_t) % do for all tdr impedances
                    param.RL_sel=port_sel(ipsl); % this used in get_TDR
                    %                     OP.interp_sparam_phase='interp_to_DC'; % better for return loss
                    %                     OP.interp_sparam_mag='trend_to_DC';
                    OP.interp_sparam_mag='linear_trend_to_DC';
                    % OP.interp_sparam_mag='extrap_to_DC_or_zero';
                    OP.interp_sparam_phase='extrap_cubic_to_dc_linear_to_inf';
                    TDR_results(izt,ipsl) = get_TDR(S, OP, param,param.Z_t(izt),ipsl);
                    if ipsl ==1
                        chdata(i).TDR11(izt).ZSR=[TDR_results(izt,1).tdr];
                        chdata(i).TDR11(izt).t=TDR_results(izt,1).t;
                        chdata(i).TDR11(izt).avgZport=[TDR_results(izt,1).avgZport];
                        if OP.PTDR, chdata(i).PDTR11(izt).ptdr=TDR_results(izt,1).ptdr_RL;end
                    else
                        chdata(i).TDR22(izt).ZSR=[TDR_results(izt,2).tdr];
                        chdata(i).TDR22(izt).t=TDR_results(izt,2).t;
                        chdata(i).TDR22(izt).avgZport=[TDR_results(izt,2).avgZport];
                        if OP.PTDR, chdata(i).PDTR22(izt).ptdr=TDR_results(izt,2).ptdr_RL;end
                    end
                    if OP.PTDR && i==1
                        if ipsl ==1
                            chdata(i).TDR11(izt).ERL=[TDR_results(izt,1).ERL];
                            chdata(i).TDR11(izt).ERLRMS=[TDR_results(izt,1).ERLRMS];
                        else
                            if ~param.FLAG.S2P
                                chdata(i).TDR22(izt).ERL=[TDR_results(izt,2).ERL];
                                chdata(i).TDR22(izt).ERLRMS=[TDR_results(izt,2).ERLRMS];
                            else
                                chdata(i).TDR22(izt).ERL=[];
                                chdata(i).TDR22(izt).ERLRMS=[];
                            end
                        end
                    else
                        chdata(i).TDR11(izt).ERL=[];
                        chdata(i).TDR22(izt).ERL=[];
                        chdata(i).TDR11(izt).ERLRMS=[];
                        chdata(i).TDR22(izt).ERLRMS=[];
                    end
                end
            end
        end
        
    end
    if OP.DISPLAY_WINDOW && OP.DEBUG && OP.TDR
        h=figure(180);set(gcf,'Tag','COM');
        if param.package_testcase_i==1 && i == 1
            if i==1
                htabgroup = uitabgroup(h);
                htab1 = uitab(htabgroup, 'Title', 'TDR TX');
                htab3 = uitab(htabgroup, 'Title', 'PTDR TX');
                hax1 = axes('Parent', htab1);
                hax3 = axes('Parent', htab3);
                if ~param.FLAG.S2P
                    htab2 = uitab(htabgroup, 'Title', 'TDR RX');
                    htab4 = uitab(htabgroup, 'Title', 'PTDR RX');
                    hax2 = axes('Parent', htab2);
                    hax4 = axes('Parent', htab4);
                end
            end
            set(h,'CurrentAxes',hax1)
            hold on
            plot(chdata(i).TDR11(izt).t(:),chdata(i).TDR11(izt).ZSR,'disp',[ chdata(i).base ' Tx port']);
            hold off
            legend (hax1, 'off');grid on;zoom xon;
            set(legend (hax1, 'show'), 'interp', 'none');
            
            if ~param.FLAG.S2P
                set(h,'CurrentAxes',hax2)
                hold on
                plot(chdata(i).TDR22(izt).t(:),chdata(i).TDR22(izt).ZSR,'disp',[ chdata(i).base ' Rx port']);
                hold off
                legend (hax2, 'off');grid on;zoom xon;
                set(legend (hax2, 'show'), 'interp', 'none');
            end
            
            set(h,'CurrentAxes',hax3)
            hold on
            if OP.PTDR
                for izt=1:length(param.Z_t)
                    msg=[ chdata(i).base ' Tx port Zt='  num2str(param.Z_t(izt),3) '  ERL=', num2str(chdata(i).TDR11(izt).ERL, 3) 'db'];
                    plot(chdata(i).TDR11(izt).t/param.ui,chdata(i).PDTR11(izt).ptdr,'disp',msg);
                    msg=['PTDR Zt=' num2str(param.Z_t(izt)/param.ui,3) ' worst sampled noise cursors ' 'Tx port Zt='  num2str(param.Z_t(izt),3) ];
                    stem(TDR_results(izt,1).WC_ptdr_samples_t/param.ui,TDR_results(izt,1).WC_ptdr_samples,'disp',msg);
                end
            end
            hold off
            legend (hax3, 'off');grid on;zoom xon;
            set(legend (hax3, 'show'), 'interp', 'none');
            if ~param.FLAG.S2P
                set(h,'CurrentAxes',hax4)
                hold on
                if OP.PTDR
                    for izt=1:length(param.Z_t)
                        msg=[ chdata(i).base ' Tx port Zt='  num2str(param.Z_t(izt),3) '  ERL=', num2str(chdata(i).TDR22(izt).ERL, 3) 'db'];
                        plot(chdata(i).TDR22(izt).t/param.ui,chdata(i).PDTR22(izt).ptdr,'disp',msg);
                        msg=['PTDR Zt=' num2str(param.Z_t(izt),3)/param.ui ' worst sampled noise cursors ' 'Rx port Zt='  num2str(param.Z_t(izt),3) ];
                        stem(TDR_results(izt,2).WC_ptdr_samples_t/param.ui,TDR_results(izt,2).WC_ptdr_samples,'disp',msg);
                    end
                end
                hold off
                legend (hax4, 'off');grid on;zoom xon;
                set(legend (hax4, 'show'), 'interp', 'none');
            end
        end
    end
    if param.FLAG.S2P, return; end
end
function S =r_parrelell2(zref,f,rpad)
S.Parameters(1,1,:) =  -zref/(rpad*(zref/rpad + 2)).*ones(1,length(f));
S.Parameters(2,2,:) =  -zref/(rpad*(zref/rpad + 2)).*ones(1,length(f));
S.Parameters(2,1,:) =   2/(zref/rpad + 2).*ones(1,length(f));
S.Parameters(1,2,:) =  2/(zref/rpad + 2).*ones(1,length(f));
% Sm=sparameters(S.Parameters,f,zref);





function [sch,schFreqAxis]=read_Nport_touchstone(touchstone_file,port_order)

%touchstone_file:  .sNp touchstone file to read
%port_order:  port reorder vector
%
%sch:  sparameter matrix
%schFreqAxis:  frequency axis

[file_path,root_name,extension]=fileparts(touchstone_file);
fid=fopen(touchstone_file);

%fetch number of ports from extension
num_ports=str2num(char(regexp(extension,'\d*','match')));

%Get option line
[optstr,opt_pos] = textscan(fid,'%s',1,'Delimiter','','CommentStyle','!');
optcell=textscan(optstr{1}{1},'%s');
optcell=optcell{1};
while isempty(optcell) || isempty(strfind(optcell{1},'#'))
    %Some touchstone files need this.  can't remember why now.  maybe lines
    %with whitespace but not empty but not commented
    [optstr,opt_pos] = textscan(fid,'%s',1,'Delimiter','','CommentStyle','!');
    optcell=textscan(optstr{1}{1},'%s');
    optcell=optcell{1};
end

%read the entire file
raw_read_data = textscan(fid,'%f %f %f %f %f %f %f %f %f','CollectOutput',true,'CommentStyle','!');
raw_column_data=raw_read_data{1};
fclose(fid);

%number of columns for 2D matrix
columns=num_ports*num_ports*2+1;

%find the frequency lines by searching for the right number of NaN
a=sum(isnan(raw_column_data),2);
if num_ports==3
    b=find(a==2);
elseif num_ports==1
    b=find(a==6);
else
    b=find(a==0);
end

num_freq=length(b);

%toss out the NaN and reshape into a 2D matrix
raw_input = raw_column_data.';
raw_input = raw_input(~isnan(raw_input));
raw_input = reshape(raw_input,columns,num_freq).';

%get the frequency mult
frequency_mult_text=optcell{2};
if(strcmpi(frequency_mult_text,'hz'))
    frequency_mult=1;
elseif(strcmpi(frequency_mult_text,'khz'))
    frequency_mult=1e3;
elseif(strcmpi(frequency_mult_text,'mhz'))
    frequency_mult=1e6;
elseif(strcmpi(frequency_mult_text,'ghz'))
    frequency_mult=1e9;
else
    error('Unsupported format for frequency multiplier %s',frequency_mult_text);
end

%get the RI/MA/DB format
format=optcell{4};
%get Z0
port_impedance=str2double(optcell(6:end))';


%grab frequency
raw_input(:,1)=raw_input(:,1)*frequency_mult;
Spar.F=raw_input(:,1);
Spar.F=transpose(Spar.F(:));


%transform data to real imaginary
%for 2.0 support, keep it in 2D form instead of 3D because we may need to process upper/lower sparse matrix definitions
if(strcmpi(format,'ri'))
    ri_data_2D=raw_input(:,2:2:end)+raw_input(:,3:2:end)*1i;
elseif(strcmpi(format,'ma'))
    mag_data=raw_input(:,2:2:end);
    rad_data=raw_input(:,3:2:end)*pi/180;
    ri_data_2D=mag_data.*cos(rad_data)+mag_data.*sin(rad_data)*1i;
elseif(strcmpi(format,'db'))
    mag_data=10.^(raw_input(:,2:2:end)/20);
    rad_data=raw_input(:,3:2:end)*pi/180;
    ri_data_2D=mag_data.*cos(rad_data)+mag_data.*sin(rad_data)*1i;
else
    error('Format %s is not supported.  Use RI MA or DB',format);
end



%transform to 3D
%allow for upper/lower matrix specification for touchstone 2.0 support
matrix_format=0;
if(matrix_format==0)
    %full
    for j=1:num_ports
        pre_out.sp(j,1:num_ports,:)=transpose(ri_data_2D( :,(j-1)*num_ports+1:j*num_ports));
    end
elseif(matrix_format==1)
    %upper
    used_ports=0;
    for j=1:num_ports
        stated_ports=num_ports-j+1;
        pre_out.sp(j,j:num_ports,:)=transpose(ri_data_2D( :,used_ports+1:used_ports+stated_ports));
        pre_out.sp(j:num_ports,j,:)=transpose(ri_data_2D( :,used_ports+1:used_ports+stated_ports));
        used_ports=used_ports+stated_ports;
    end
elseif(matrix_format==2)
    %lower
    used_ports=0;
    for j=1:num_ports
        stated_ports=j;
        pre_out.sp(j,1:j,:)=transpose(ri_data_2D( :,used_ports+1:used_ports+stated_ports));
        pre_out.sp(1:j,j,:)=transpose(ri_data_2D( :,used_ports+1:used_ports+stated_ports));
        used_ports=used_ports+stated_ports;
    end
else
    error('Matrix format is not supported.  Use Full, Lower, or Upper');
end


%check for swapping the 2 port matrix (required on 1.x spec)
two_port_swap=1;
if(num_ports==2 && two_port_swap==1)
    temp=pre_out.sp(1,2,:);
    pre_out.sp(1,2,:)=pre_out.sp(2,1,:);
    pre_out.sp(2,1,:)=temp;
end

Spar.S=pre_out.sp;
Spar.Z0=transpose(port_impedance(:));

if length(Spar.Z0)>1
    error('Each port must have the same reference impedance');
end
if ~isequal(Spar.Z0,50)
    warning('Reference impedance of %0.6g ohms renormalized to 50 ohms',Spar.Z0);
    %Renormalize to 50 ohms
    rho=(50-Spar.Z0)/(50+Spar.Z0);
    p=num_ports;
    s_old=Spar.S;
    for k=1:num_freq
        Spar.S(:,:,k)=inv(eye(p,p)-rho*s_old(:,:,k))*(s_old(:,:,k)-rho*eye(p,p));
    end
end

%These operations sync up with COM style Spar matrix
%1:  put frequency as first dimension
sch=shiftdim(Spar.S,2);
%2:  reorder ports according to "ports" input
sch=sch(:,port_order,port_order);
schFreqAxis=Spar.F;
function [chdata, param] = read_PR_files(param, OP, chdata)
%% Read in pulse response
% extract s-parameter data from files and apply tx and rx filters as well as package filters
num_files=length(chdata);
M=param.samples_per_ui;
if ~OP.DISPLAY_WINDOW, fprintf('reading file '); end
for i=1:num_files
    if OP.DISPLAY_WINDOW; hwaitbar=waitbar(0);end
    progress = i/num_files;
    if OP.DISPLAY_WINDOW
        [~,a]=fileparts(chdata(i).filename);
        waitbar(progress, hwaitbar, ['Processing ' a]); figure(hwaitbar); drawnow;
    else
        fprintf('%i ',i);
    end
    switch chdata(i).ext
        case '.csv'
            vt=load(chdata(i).filename);
            % many PR's need to have precursors added: we'll a 2 here RIM 5-24-2021
            chdata(i).uneq_pulse_response=[ zeros(3*param.samples_per_ui, 1) ; vt(:,2) ];
            dt=vt(2,1)-vt(1,1);
            chdata(i).t= [ (0:3*param.samples_per_ui-1).'*dt ; vt(:,1)+3*param.samples_per_ui*dt];
             
            
            step_shifting_vector=kron(ones(1,floor(length( chdata(i).uneq_pulse_response)/M)) ,[ 1  zeros(1,M-1) ]) ;
            step_response=filter(step_shifting_vector,1,chdata(i).uneq_pulse_response);
            Vf=step_response(end);
            chdata(i).uneq_imp_response=step_response-circshift(step_response,1); % too noisey to be usefull
            chdata(i).uneq_imp_response(1)=chdata(i).uneq_imp_response(2);
            
    end
end
function [param,OP]= read_ParamConfigFile(paramFile,OP)
%warning('off','MATLAB:xlsread:Mode'); % suppress warning messages for reading the settings file from XLS
[filepath,name,ext] = fileparts(paramFile);
if ~isempty(filepath)
    filepath=[filepath '\'];
end
matcongfile=[filepath name '.mat'];
try
    switch upper(ext)
        case upper('.mat')
            load(matcongfile)
        case upper('.csv')
            [na1, na2, parameter] = xlsread(paramFile);
        otherwise
            [na1, na2, parameter] = xlsread(paramFile,'COM_Settings','',''); %#ok<ASGLU> % Import data from the settings file (imports the entire sheet)
    end
    
catch ME %#ok<NASGU>
    warning('off','MATLAB:xlsread:Mode'); % suppress warning messages for reading the settings file from XLS
    switch upper(ext)
        case upper('.mat')
            load(matcongfile)
        case upper('.csv')
            [na1, na2, parameter] = xlsread(paramFile);
        otherwise
            [na1, na2, parameter] = xlsread(paramFile,'COM_Settings','','basic'); %#ok<ASGLU> % Import data from the settings file (imports the entire sheet)
    end
end

%% New section to parse .START package data
first_column_data = parameter(:,1);
start_data_rows = find(strcmp(first_column_data,'.START'));
if ~isempty(start_data_rows)
    end_data_rows = find(strcmp(first_column_data,'.END'));
    if length(start_data_rows) ~= length(end_data_rows)
        error('Number of .START and .END must be the same');
    end
    first_start_row = start_data_rows(1);
    special_parameter = parameter;
    parameter = parameter(1:first_start_row-1,:);
    for j=1:length(start_data_rows)
        this_block = special_parameter(start_data_rows(j)+1:end_data_rows(j)-1,:);
        pkg_name = special_parameter{start_data_rows(j),2};
        
        %Read all the parameters that make up a package
        PKG_param = read_package_parameters(this_block);
        
        %save the data in a field revealed by pkg_name
        param.PKG.(pkg_name) = PKG_param;
        
        
    end
end
%Allow specification of TX and RX package section through PKG_NAME keyword
%the values must match package blocks specified in .START sections
param.PKG_NAME= xls_parameter(parameter, 'PKG_NAME', false,'');
if isnan(param.PKG_NAME)
    param.PKG_NAME = '';
end
if isempty(param.PKG_NAME)
    param.PKG_NAME = {};
else
    param.PKG_NAME = strsplit(param.PKG_NAME);
end
if ~isempty(param.PKG_NAME) && ~isfield(param,'PKG')
    error('PKG_NAME can only be used if .START blocks for package parameters are used');
end
for j=1:length(param.PKG_NAME)
    if ~isfield(param.PKG,param.PKG_NAME{j})
        error('Package Block "%s" not found',param.PKG_NAME{j});
    end
end

%%
% just need to define so we can pass
param.c=[.4e-12 .4e-12]; 
param.alen=[ 20 30 550  ];
param.az=[100 120 100];

% make control for package/channel reflection control
param.kappa1= xls_parameter(parameter, 'kappa1', true,1); % if set 0 reflection at tp0 are omitted from COM
param.kappa2= xls_parameter(parameter, 'kappa2', true,1); % if set 0 reflection at tp5 are omitted from COM


% make compatible with presentation of kappa

% Default values are given for parameters when they are common to all clauses in 802.3bj and 803.2bm.

OP.dynamic_txffe = xls_parameter(parameter, 'Dynamic TXFFE', false,1); % Temporary switch for testing new optimize_fom with dynamic txffe
OP.FloatingDFE_Development = xls_parameter(parameter, 'FloatingDFE_Development', false,1); % Temporary switch for testing new floating dfe routine

param.fb = xls_parameter(parameter, 'f_b')*1e9; % Baud (Signaling) rate in Gbaud
param.f2 = xls_parameter(parameter, 'f_2', true, param.fb/1e9 )*1e9; % frequency in GHz for intergration compuation of ICN or FOM_Ild in GHz
param.max_start_freq = xls_parameter(parameter, 'f_min')*1e9; % minimum required frequency start for s parameters
param.f1 = xls_parameter(parameter, 'f_1', true, param.max_start_freq/1e9)*1e9; % start frequency for ICN and ILD calculations in GHz
param.max_freq_step = xls_parameter(parameter, 'Delta_f')*1e9; % freqency step
param.tx_ffe_c0_min = xls_parameter(parameter, 'c(0)', false); % TX equalizer cursor minimum value (actual value is calculated as 1-sum(abs(tap)), Grid seat ingored for when C(0) is below this value

if OP.dynamic_txffe
    found_pre=1;
    pre_count=1;
    while found_pre
        [p,found_pre]=xls_parameter_txffe(parameter,sprintf('c(-%d)',pre_count));
        if found_pre
            field_name=sprintf('tx_ffe_cm%d_values',pre_count);
            param.(field_name)=p;
            pre_count=pre_count+1;
        end
    end
    found_post=1;
    post_count=1;
    while found_post
        [p,found_post]=xls_parameter_txffe(parameter,sprintf('c(%d)',post_count));
        if found_post
            field_name=sprintf('tx_ffe_cp%d_values',post_count);
            param.(field_name)=p;
            post_count=post_count+1;
        end
    end
else
param.tx_ffe_cm4_values = xls_parameter(parameter, 'c(-4)', true,0); % TX equalizer pre cursor tap -4 individual settings or range. If not present ignored
param.tx_ffe_cm3_values = xls_parameter(parameter, 'c(-3)', true,0); % TX equalizer pre cursor tap -3 individual settings or range. If not present ignored
param.tx_ffe_cm2_values = xls_parameter(parameter, 'c(-2)', true,0); % TX equalizer pre cursor tap -2 individual settings or range. If not present ignored
param.tx_ffe_cm1_values = xls_parameter(parameter, 'c(-1)', true,0); % TX equalizer pre cursor tap -1 individual settings or range. If not present ignored
param.tx_ffe_cp1_values = xls_parameter(parameter, 'c(1)', true,0); % TX equalizer post cursor tap 1  individual settings or range. If not present ignored
param.tx_ffe_cp2_values = xls_parameter(parameter, 'c(2)', true,0); % TX equalizer post cursor tap 2  individual settings or range. If not present ignored
param.tx_ffe_cp3_values = xls_parameter(parameter, 'c(3)', true,0); % TX equalizer post cursor tap 3  individual settings or range. If not present ignored
end
param.ndfe = xls_parameter(parameter, 'N_b'); % Decision feedback fixed equalizer (DFE) length
param.N_v = xls_parameter(parameter, 'N_v',true,param.ndfe); % number of UI used to compute Vf
param.D_p = xls_parameter(parameter, 'D_p',true, 4 ); % number of precursor UI's used to compute Vf Default to 10
param.N_bx = xls_parameter(parameter, 'N_bx', true, param.ndfe ); % Used for ERL to Compensate for a number of Ui assoicated with the DFE
% support for floating taps
param.N_bg=xls_parameter(parameter, 'N_bg', true,0); % number of group of floating tap. Used as a switch, 0 means no float
param.N_bf=xls_parameter(parameter, 'N_bf', true,6); % number of taps in group
param.N_bmax=xls_parameter(parameter, 'N_bmax', true, param.ndfe); % UI span for floating taps
param.N_bmax=xls_parameter(parameter, 'N_f', true, param.ndfe); % UI span for floating taps. replaced by N_bmax
param.N_f=xls_parameter(parameter, 'N_f', true, param.ndfe); % UI span for floating taps for rX FFE
if param.N_bg == 0, param.N_bmax=param.ndfe; end
param.bmaxg=xls_parameter(parameter, 'bmaxg', true,0.2); % max DFE value for floating taps

% support for tail tap power limitations
param.B_float_RSS_MAX=xls_parameter(parameter, 'B_float_RSS_MAX', true,0); % floating DFE tap start for RSS floating tap limit
param.N_tail_start=xls_parameter(parameter, 'N_tail_start', true,0); % start range for max RSS limit for DFE taps
%
param.dfe_delta = xls_parameter(parameter, 'N_b_step', true,0); % discreatiztion of DFE, 0 disables and is not normally used
param.ffe_pre_tap_len=xls_parameter(parameter, 'ffe_pre_tap_len', true,0); % RX ffe pre cursor tap length
param.RxFFE_cmx=param.ffe_pre_tap_len;
param.ffe_post_tap_len=xls_parameter(parameter, 'ffe_post_tap_len', true,0); % Rx FFE post cursor tap length
param.RxFFE_cpx=param.ffe_post_tap_len;
param.ffe_tap_step_size=xls_parameter(parameter, 'ffe_tap_step_size', true,0); % Rx FFE tap step size
param.RxFFE_stepz=param.ffe_tap_step_size;
param.ffe_main_cursor_min=xls_parameter(parameter, 'ffe_main_cursor_min', true,1); % Rx FFE main cursor miminum 
param.ffe_pre_tap1_max=xls_parameter(parameter, 'ffe_pre_tap1_max', true,.7); % Rx FFE precursor tap1 limit
param.ffe_post_tap1_max=xls_parameter(parameter, 'ffe_post_tap1_max', true,.7); % Rx FFE post cursor tap1 limit
param.ffe_tapn_max=xls_parameter(parameter, 'ffe_tapn_max', true,.7); % Rx FFE precursor tapn limit
param.ffe_backoff=xls_parameter(parameter, 'ffe_backoff', true,0); % see if better zero foreced solution is better by backing off the number specified FFE taps one at at time
if param.RxFFE_cmx ~= 0 || param.RxFFE_cpx ~=0
    OP.RxFFE= true;
else
    OP.RxFFE=false;
end
param.num_ui_RXFF_noise=xls_parameter(parameter, 'num_ui_RXFF_noise', true,2048); % Rx FFE precursor tapn limit

param.g_DC_HP_values = xls_parameter(parameter, 'g_DC_HP', true,[]); % CTF AC-DC gain list (GDC2)
param.f_HP = 1e9*xls_parameter(parameter, 'f_HP_PZ', true, []); % CFT pole pole zero pair in GHz for low frequency CTF
param.f_HP_Z = 1e9*xls_parameter(parameter, 'f_HP_Z', true, []); % CFT zero fz1 is in GHz. Normally a list for 120e. Not normally use elsewise
param.f_HP_P = 1e9*xls_parameter(parameter, 'f_HP_P', true, []); % CFT pole fp2 is in GHz. Normally a list for 120e. Not normally use elsewise


param.Min_VEO= xls_parameter(parameter, 'EH_min', true,0); % used when PMD_type is C2M
param.Max_VEO= xls_parameter(parameter, 'EH_max', true,inf); % used when PMD_type is C2M and is not really computed per spec
param.Min_VEO_Test= xls_parameter(parameter, 'EH_min_test', true,0); % Older syntax. Used when PMD_type is C2M. This allow EH to go below EH_min. If set to zero it is ignored (same as Min_VEO_test)
param.Min_VEO_Test= xls_parameter(parameter, 'Min_VEO_Test', true,param.Min_VEO_Test); % used when PMD_type is C2M. This allow EH to go blow EH_min. If set to Zero it is ignored

param.CTLE_type= xls_parameter(parameter, 'CTLE_type', false,'CL93'); % Sets the CTLE type default is poles and zeros (i.e. not a list of poles as in 120e) 
if ~isempty(param.g_DC_HP_values) ; param.CTLE_type='CL120d';end % overrides CL93 if g_DC_HP_values are a spreadsheet entry. Mostly used when baud rare is >= 50Gbps
if ~isempty(param.f_HP_Z) ; param.CTLE_type='CL120e';end % overrides CL93 and CL120d if f_HP_Z is a spreadsheet entry
% always read in main ctle values. They would be interpreted different baseed
% on the clause they apply because of different CTF equations
param.ctle_gdc_values = xls_parameter(parameter, 'g_DC', true); % AC-DC gain list
param.CTLE_fp1 = 1e9*xls_parameter(parameter, 'f_p1', true, param.fb/4); % fp1 is in GHz
param.CTLE_fp2 = 1e9*xls_parameter(parameter, 'f_p2', true, param.fb); % fp2 is in GHz
param.CTLE_fz = 1e9*xls_parameter(parameter, 'f_z', true, param.fb/4); % fz is in GHz
% the contex of the poles an zeros are determined by the clause
switch param.CTLE_type
    case 'CL93'
        param.ctle_gdc_values = xls_parameter(parameter, 'g_DC', true); % Continuous time filter DC gain settings (G_DC) or range as specified in Annex 93A
        param.CTLE_fp1 = 1e9*xls_parameter(parameter, 'f_p1', true, param.fb/4); % fp1 is in GHz
        param.CTLE_fp2 = 1e9*xls_parameter(parameter, 'f_p2', true, param.fb); % fp2 is in GHz
        param.CTLE_fz = 1e9*xls_parameter(parameter, 'f_z', true, param.fb/4); % fz is in GHz
    case 'CL120d'
        param.g_DC_HP_values = xls_parameter(parameter, 'g_DC_HP', true,[]); % Continuous time filter DC gain settings (G_DC2)
        param.f_HP = 1e9*xls_parameter(parameter, 'f_HP_PZ', true, []); % fLF is in GHz
    case 'CL120e'
        % re adjust to get TD_CTLE to work with C:120e equation without
        % changing TD_CTLE code
        param.CTLE_fz =param.CTLE_fz ./ 10.^(param.ctle_gdc_values/20);
end
param.GDC_MIN = xls_parameter(parameter, 'GDC_MIN',true, 0); % max ACDC gain, if 0 ignore
param.cursor_gain=xls_parameter(parameter, 'crusor_gain', true,0); % only FFE and not supported
param.a_thru = xls_parameter(parameter, 'A_v', true); % Victim differential peak source output voltage (half of peak to peak)
param.a_fext = xls_parameter(parameter, 'A_fe', true); % FEXT aggressor differential peak source output voltage (half of peak to peak)
param.a_next = xls_parameter(parameter, 'A_ne', true); % NEXT aggressor differential peak source output voltage (half of peak to peak)
param.a_icn_fext = xls_parameter(parameter, 'A_ft', true, param.a_fext); % FEXT aggressor amplitude for ICN. Defaults to A_fe if not specified
param.a_icn_next = xls_parameter(parameter, 'A_nt', true, param.a_next );% NEXT aggressor amplitude for ICN. Defaults to A_ne if not specified
param.levels = xls_parameter(parameter, 'L'); % number of symbols levels (PAM-4 is 4, NRZ is 2)
param.specBER = xls_parameter(parameter, 'DER_0'); % Target detector error ratio
param.pass_threshold = xls_parameter(parameter, 'COM Pass threshold',false,0); % the pass fail threshold for COM in dB
param.ERL_pass_threshold = xls_parameter(parameter, 'ERL Pass threshold',false,0); % the pass fail threshold for ERL in dB
param.VEC_pass_threshold = xls_parameter(parameter, 'VEC Pass threshold',false,0);% the pass fail threshold for VEC in dB only used when PMD_type is C2M

param.sigma_RJ = xls_parameter(parameter, 'sigma_RJ'); % rms of of random jitter
param.A_DD = xls_parameter(parameter, 'A_DD'); % Normalized peak dual-Dirac noise, this is half of the total bound uncorrelated jitter (BUJ) in UI
param.eta_0 = xls_parameter(parameter, 'eta_0'); % One-sided noise spectral density (V^2/GHz).Input refered noise at TP5. Input referred noise at TP5
param.SNDR = xls_parameter(parameter, 'SNR_TX', true); % Transmitter SNDR noise in dB
param.R_LM = xls_parameter(parameter, 'R_LM'); % Ratio of level separation mismatch. Relevant when not PAM-2 (NRZ).
param.samples_per_ui = xls_parameter(parameter, 'M', 32); % Samples per UI
param.ts_sample_adj_range = xls_parameter(parameter, 'sample_adjustment', true, [0 0]); %sample point adjustment range
param.ts_anchor = xls_parameter(parameter, 'ts_anchor', true, 0); %choice of sampling routine.  0=MM, 1=Peak, 2=max DV (max cursor minus precursor)
% This will keep bmax length 0 if Nb=0

%AJG021820
param.bmax(1:param.ndfe)     = xls_parameter(parameter, 'b_max(1)'); % DFE magnitude limit, first coefficient(ignored if Nb=0)
if isempty(param.bmax)
    param.bmin=param.bmax;
else
    param.bmin(1:param.ndfe)     =xls_parameter(parameter, 'b_min(1)', true,-param.bmax(1) ); % DFE negative magnitude limit. If not specified it defaults to -bmax.

end
if param.ndfe >= 2
    param.bmax(2:param.ndfe) = xls_parameter(parameter, 'b_max(2..N_b)', true, .2); % DFE magnitude limit, second coefficient and on (ignored if Nb<2). Can be a regualar expression
    param.bmin(2:param.ndfe) = xls_parameter(parameter, 'b_min(2..N_b)', true, -1*param.bmax(2:param.ndfe) ); % DFE negative magnitude limit, if not specified it defaults to -b_max(2..N_b)
end

param.gqual=xls_parameter(parameter, 'G_Qual', true,[]);% G_Qual are the dB ranges of g_DC g DC )which correspond tog_DC_HP (g DC2)
param.g2qual=xls_parameter(parameter, 'G2_Qual', true,[]); % G2_Qual limit values of g_DC_HP (g DC2 ) which corresponds to ranges of g_DC g DC specified with G_QUAL
%verify gqual and gqual2 input
if ~isempty(param.gqual) || ~isempty(param.g2qual)
    if size(param.gqual,1)~=length(param.g2qual)
        error('gqual and g2qual size mismatch');
    end
    if size(param.gqual,2)~=2
        error('gqual must be Nx2 matrix');
    end
end


% eval if string for all three - can use different for TX and RX
param.C_pkg_board = xls_parameter(parameter, 'C_p', true)*1e-9; % C_p in nF (single sided)
param.C_diepad = xls_parameter(parameter, 'C_d', true)*1e-9; % C_d in nF (single sided)
% [ahealey] Read values for optional compensating L and "bump" C
param.L_comp = xls_parameter(parameter, 'L_s', true, 0)*1e-9; % L_s in nH (single sided)
param.C_bump = xls_parameter(parameter, 'C_b', true, 0)*1e-9; % C_b in nF (single sided)
% [ahealey] End of modifications.
param.C_v = xls_parameter(parameter, 'C_v', true,0)*1e-9; % C_v in nF (via cap)  (single sided)
param.R_diepad = xls_parameter(parameter, 'R_d', true); % Die source termination resistance  (single sided)
param.Z_t = xls_parameter(parameter, 'Z_t', true,50); %  single sided source termination reference resistance for TDR and ERL 
param.TR_TDR = xls_parameter(parameter, 'TR_TDR', true , 8e-3); %  Gaussian shaped transition time for TDR source in ns 


param.Z0 = xls_parameter(parameter, 'R_0', 50); % 
param.z_p_tx_cases = xls_parameter(parameter, 'z_p (TX)', true).'; % List of victim transmitter package trace lengths in mm, one per case
[ncases, mele]=size(param.z_p_tx_cases);
if mele ==2
    param.flex=2;
elseif mele==4
    param.flex=4;
elseif mele==1
    param.flex=1;
else
    error(sprintf('config file syntax error'))
end

% board parameters
param.C_0 = xls_parameter(parameter, 'C_0', true,0)*1e-9; % If Include PCB is set to 1, near device single ended capacitance C0  in nF is added  
param.C_1 = xls_parameter(parameter, 'C_1', true,0)*1e-9; % if Include PCB is set to 1, connector side single ended capacitance C1 in nF is added 
%
param.z_p_next_cases = xls_parameter(parameter, 'z_p (NEXT)', true).'; % List of NEXT transmitter package trace lengths in mm, one per case
[ncases1, mele1]=size(param.z_p_next_cases);
if ncases ~= ncases1 || mele ~= mele1
    error('All TX, NEXT, FEXT, Rx cases must agree');
else
end
param.z_p_fext_cases = xls_parameter(parameter, 'z_p (FEXT)', true).'; % List of FEXT transmitter package trace lengths in mm, one per case
[ncases1, mele1]=size(param.z_p_fext_cases);
if ncases ~= ncases1 ||  mele ~= mele1
    error('All TX, NEXT, FEXT, Rx cases must agree');
else
end
param.z_p_rx_cases = xls_parameter(parameter, 'z_p (RX)', true).'; % List of FEXT receiver package trace lengths in mm, one per case
[ncases1, mele1]=size(param.z_p_rx_cases);
if ncases ~= ncases1 ||  mele ~= mele1
    error('All TX, NEXT, FEXT, Rx cases must agree');
else
end
% Table 93A-3 parameters
param.pkg_gamma0_a1_a2 = xls_parameter(parameter, 'package_tl_gamma0_a1_a2', true, [0 1.734e-3 1.455e-4]); %Fitting parameters for package model per unit length. First element is in 1/mm and affects DC loss of package model . Second element is in ns1/2/mm and affects loss proportional to sqrt(f). Third element is in ns/mm and affects loss proportional to f.
param.pkg_tau = xls_parameter(parameter, 'package_tl_tau', true, 6.141e-3); % Package model transmission line delay ns/mm
param.pkg_Z_c = xls_parameter(parameter, 'package_Z_c', true, 78.2).';% Package model transmission line characteristic impedance [ Tx , Rx ]
[ ncases1, mele1]=size(param.pkg_Z_c);% 
if   mele ~= mele1
    error('tx rx pairs must have thesame number element entries as TX, NEXT, FEXT, Rx');
else
end
if mele1==2 % fuill in a array if only a 2 element flex package is specified
    for ii=1:ncases
        param.z_p_fext_casesx(ii,:)=  [param.z_p_fext_cases(ii,:)' ;[ 0 ; 0 ]]';
        param.z_p_next_casesx(ii,:)=  [param.z_p_next_cases(ii,:)' ;[ 0 ; 0 ]]';
        param.z_p_tx_casesx(ii,:)=  [param.z_p_tx_cases(ii,:)' ;[ 0 ; 0 ]]';
        param.z_p_rx_casesx(ii,:)=  [param.z_p_rx_cases(ii,:)' ;[ 0 ; 0 ]]';
    end
    param.z_p_fext_cases =  param.z_p_fext_casesx;
    param.z_p_next_cases=  param.z_p_next_casesx;
    param.z_p_tx_cases=  param.z_p_tx_casesx;
    param.z_p_rx_cases=  param.z_p_rx_casesx;
    param.pkg_Z_c=[param.pkg_Z_c' ;[ 100 100 ; 100 100 ]]';
end
param.PKG_Tx_FFE_preset =xls_parameter(parameter, 'PKG_Tx_FFE_preset', true, 0); % RIM 08-18-2022 for Tx preset capability

% Table 92-12 parameters
param.brd_gamma0_a1_a2 = xls_parameter(parameter, 'board_tl_gamma0_a1_a2', true, [0 4.114e-4 2.547e-4]); % Fitting parameters for package model per unit length. First element is in 1/mm and affects DC loss of package model . Second element is in ns1/2/mm and affects loss proportional to sqrt(f). Third element is in ns/mm and affects loss proportional to f.
param.brd_tau = xls_parameter(parameter, 'board_tl_tau', true, 6.191e-3);% Board model transmission line delay ns/mm
param.brd_Z_c = xls_parameter(parameter, 'board_Z_c', true, 109.8); % Board model transmission line characteristic impedance [ Tx , Rx ]
param.z_bp_tx = xls_parameter(parameter, 'z_bp (TX)', true, 151); %  Victim transmitter board trace lengths in mm
param.z_bp_next = xls_parameter(parameter, 'z_bp (NEXT)', true, 72);% Next Assessor transmitter board trace lengths in mm
param.z_bp_fext = xls_parameter(parameter, 'z_bp (FEXT)', true, 72);% Rext Assessor transmitter board trace lengths in mm
param.z_bp_rx = xls_parameter(parameter, 'z_bp (RX)', true, 151);% Victim receiver board trace lengths in mm

% Unofficial parameters
param.snpPortsOrder = xls_parameter(parameter, 'Port Order', true, [1 3 2 4]); % s parameter port order [ tx+ tx- rx+ rx-]
param.delta_IL=xls_parameter(parameter, 'delta_IL', false, 1); % experiemnal
% Deprecated parameters - affect only frequency domain analysis.
param.f_v = xls_parameter(parameter, 'f_v', true, 4); % For FOM_ILD: Transiton rate cut off frequency for ICN/ILD calc in terms of fb
param.f_f = xls_parameter(parameter, 'f_f', true, 4); % For ICN: Fext transiton rate cut off frequency for ICN calc in terms of fb
param.f_n = xls_parameter(parameter, 'f_n', true, 4); % For ICN: Next transiton rate cut off frequency for ICN calc in terms of fb
param.f_r = xls_parameter(parameter, 'f_r', true, 4); % reference receive filter in COM and in ICN/FOM_ILD calcs in terms of fb
param.fb_BT_cutoff= xls_parameter(parameter, 'TDR_f_BT_3db', true, 0.4730); % Bessel-Thomson 3 dB cut off freqeuncy in terms of fb
param.BTorder =  xls_parameter(parameter, 'BTorder', false, 4); % Bessel function order
param.RC_Start =  xls_parameter(parameter, 'RC_Start', false, param.fb/2); % start frequency for raised cosine filter
param.RC_end =  xls_parameter(parameter, 'RC_end', false, param.fb*param.f_r ); % end frequency for raised cosine filter
param.beta_x= xls_parameter(parameter, 'beta_x', false, 0);%  (for ERL) use default
param.rho_x= xls_parameter(parameter, 'rho_x', false, .618); % (for ERL) use default
param.tfx= xls_parameter(parameter, 'fixture delay time', true, -1);% fixture delay time (for ERL)
param.Grr_limit=xls_parameter(parameter, 'Grr_limit', false, 1); % either do no use or set to 1 (for ERL)
param.Grr=xls_parameter(parameter, 'Grr', false, param.Grr_limit);% either do no use or set to 1 (for ERL)
param.Gx=xls_parameter(parameter, 'Gx', false, 0); % ERL parameter param.Grr, This is used is the COM code
switch param.Gx
    case 0
        param.Grr=param.Grr; % just use older Grr ir gx not specified
    case 1
        param.Grr=2; % use newer Grr
end

param.LOCAL_SEARCH=xls_parameter(parameter,'Local Search',true,0); % Decreases COM compute time. Aetting to 2 seems ok ,if 0 search is full grid
% Operational control variables
%OP.include_pcb = xls_parameter(parameter, 'Include PCB (table 92-13)', false, 0);
param.Tukey_Window=xls_parameter(parameter,'Tukey_Window',true,0); % required for ERL. Set to 1. Default is 0.
param.Noise_Crest_Factor= xls_parameter(parameter, 'Noise_Crest_Factor', true, 0); % Normally not used. If set this is q factor used for quantized Gaussian PDFs
param.AC_CM_RMS = xls_parameter(parameter, 'AC_CM_RMS', true, 0); % AC_CM_RMS is the CM BBN AWGN RMS at COM source point. Default is 0. Adds common mode noise source to the COM signal path for the through channel
param.ACCM_MAX_Freq=xls_parameter(parameter, 'ACCM_MAX_Freq', true, param.fb); % F max for integrating ACCM voltage in Hz. Default is fb
param.T_O = xls_parameter(parameter, 'T_O', true, 0 ); % Units are mUI. Histogram for VEC and VEO are computed over T_s +/- T_O.  
param.T_O = xls_parameter(parameter, 'T_h', true, param.T_O  ); % superceded with T_O but is the internal values that is used. Do not use.

param.samples_for_C2M =xls_parameter(parameter, 'samples_for_C2M', true, 100 ); % Finer sampling in terms of samples per UI for c2m histgram analysis.

OP.Histogram_Window_Weight=xls_parameter(parameter, 'Histogram_Window_Weight', false, 'rectangle' );  %Weighting for VEC and VEO are histogram processing. Type are Gaussian,Dual Rayleigh,Triangle, and Rectangle (default)
param.sigma_r=xls_parameter(parameter, 'sigma_r', true, .020 ); % sigma_r for 0.3ck Gaussian histogram window. Unit are UI. Preferred usage.
param.Qr=xls_parameter(parameter, 'Qr', true, param.sigma_r ); % sigma_r replaces Qr gasussian histogram window. Unit are UI
param.QL=xls_parameter(parameter, 'QL', true, param.T_O/param.Qr/1000 ); % superceded with sigma_r but is the internal values that is used

%%

param.skew_ps=xls_parameter(parameter, 'skew_ps', true, 0 );% experiment p/n skew. Not used.
param.imb_Z_fctr=xls_parameter(parameter, 'imb_Z_fctr', true, 1 ); % exprimental p/n impedance missmatch.  Not used.
param.imb_C_fctr=xls_parameter(parameter, 'imb_C_fctr', true, 1 ); % exprimental p/n capacitance missmatch.  Not used.
param.awgn_mv=param.AC_CM_RMS;
param.flip=xls_parameter(parameter, 'flip', true, 0 );  % exprimental p/n missmatch flip.  Not used.
param.f_hp=xls_parameter(parameter, 'f_hp', true, 0 );  % for rx testing for eq 162-12 if 0 (default) then rx test using rx bbn 


%% Adding new parameters to reveal whether Floating DFE or Floating RXFFE is used
% This removes the dependency on checking param.N_bg (that is no longer valid to reveal if floating DFE is used)
param.Floating_RXFFE=false;
param.Floating_DFE=false;
if param.N_bg > 0
    param.Floating_DFE=true;
end
if OP.RxFFE
    param.Floating_DFE=false;
    if param.N_bg > 0
        param.Floating_RXFFE=true;
    end
end
%% for introducing Tx or Rx skew on p leg or n leg
param.Txpskew=xls_parameter(parameter, 'Txpskew', true, 0 );  % Tx p skew in ps
param.Txnskew=xls_parameter(parameter, 'Txnskew', true, 0 );  % Tx n skew in ps
param.Rxpskew=xls_parameter(parameter, 'Rxpskew', true, 0 );  % Rx p skew in ps
param.Rxnskew=xls_parameter(parameter, 'Rxnskew', true, 0 );  % Rx n skew in ps

%%
OP.include_pcb = xls_parameter(parameter, 'Include PCB', false); % Used to add a PCB one each side of the passed s-parameters.
OP.exit_if_deployed = xls_parameter(parameter, 'exit if deployed', false,0); % may need set when COM is an exe
OP.INCLUDE_CTLE = xls_parameter(parameter, 'INCLUDE_CTLE', false, 1); % do not use 
OP.EXE_MODE= xls_parameter(parameter, 'EXE_MODE', false, 1);% 12/21 0:legacy 1:fast 2:superfast default is 1.
OP.INCLUDE_FILTER = xls_parameter(parameter, 'INCLUDE_TX_RX_FILTER', false, 1); % do not use
OP.force_pdf_bin_size = xls_parameter(parameter, 'Force PDF bin size', false, 0); % do not use
OP.BinSize = xls_parameter(parameter, 'PDF bin size', false, 1e-5); % set lower for faster computation time but less accuracy. 
OP.DEBUG = xls_parameter(parameter, 'DIAGNOSTICS', false, false); % supresss some interim compuation value printouts
OP.DISPLAY_WINDOW = xls_parameter(parameter, 'DISPLAY_WINDOW', false, true); % controls if graph plots are displayed. Typically goes along with DIAGNOSTICS
OP.CSV_REPORT = xls_parameter(parameter, 'CSV_REPORT', false, true); % saves all the output parameters to a CSV file in the results directory, If DIAGNOSTICS is set then a mat file is also created
OP.SAVE_TD=xls_parameter(parameter, 'SAVE_TD', false, false); % Save the time domian waveforms. FIR, PR etc. in an output structure
OP.SAVE_FIGURES=xls_parameter(parameter, 'SAVE_FIGURES', false, false); % save displayed figures in the results directory
OP.SAVE_FIGURE_to_CSV=xls_parameter(parameter, 'SAVE_FIGURE_to_CSV', false, false); % does not work. do not use.
OP.GET_FD = xls_parameter(parameter, 'Display frequency domain', false, OP.GET_FD); % Not normally set in the config file. It is normally just set to true to get FD plots
OP.INC_PACKAGE = xls_parameter(parameter, 'INC_PACKAGE', false, true); % warning: INC_PACKAGE=0 not fully supported, instead, set Zp,Cd, and Cp parameters to zero and Zp select to 1
if ~OP.INC_PACKAGE
    fprintf('<strong> Warning!!! INC_PACKAGE=0 not fully supported, instead, set Zp,Cd, and Cp parameters to zero and Zp select to 1 </strong>\n');
end

OP.EW = xls_parameter(parameter, 'EW', false, false); % RIM 3-18-2021 change defaults
OP.IDEAL_TX_TERM = xls_parameter(parameter, 'IDEAL_TX_TERM', false, false);
if OP.IDEAL_TX_TERM
    fprintf('<strong> Warning!!! IDEAL_TX_TERM not supported, instead, set Zp,Cd, and Cp parameters to zero and Zp select to 1 </strong>\n');
end
OP.IDEAL_RX_TERM = xls_parameter(parameter, 'IDEAL_RX_TERM', false, false);
if OP.IDEAL_RX_TERM
    fprintf('<strong> Warning!!! IDEAL_RX_TERM not supported, instead, set Zp,Cd, and Cp parameters to zero and Zp select to 1 </strong>\n');
end

OP.TDMODE = xls_parameter(parameter, 'TDMODE',false, OP.TDMODE);   % Enables the the use of pulse response instead of s-parameters. Assumes no packages or the packages are included in the PR. Default is 0.

OP.FT_COOP = xls_parameter(parameter, 'FT_COOP',false, false); % obsolete do not use.
OP.RESULT_DIR = regexprep(xls_parameter(parameter, 'RESULT_DIR'), '\\', filesep); % directory where results like csv, mat, and/or figure files will be written
OP.RESULT_DIR=strrep(OP.RESULT_DIR,'{date}',date);
OP.BREAD_CRUMBS = xls_parameter(parameter, 'BREAD_CRUMBS',false, false); % if DIAGNOSTICS is set then param, OP, and chdata are include in the output for each run
OP.BREAD_CRUMBS_FIELDS = xls_parameter(parameter, 'BREAD_CRUMBS_FIELDS',false, ''); % if BREAD_CRUMBs is enabled, this file controls what chdata fields are included
OP.COM_CONTRIBUTION_CURVES = xls_parameter(parameter, 'COM_CONTRIBUTION',false,0);   % Default is 0. If set to 1 then a bar graph of COM contributors is produce instead of bathtub curves
OP.ENFORCE_CAUSALITY = xls_parameter(parameter, 'Enforce Causality', false, 0);% default is 0. Not recommended
OP.EC_REL_TOL = xls_parameter(parameter, 'Enforce Causality REL_TOL', false, 1e-2); % Relative Tolerance parameter for causality, Hard enforcement, 1e-3, Soft enforcement,  1e-2
OP.EC_DIFF_TOL = xls_parameter(parameter, 'Enforce Causality DIFF_TOL', false, 1e-3); % Difference Tolerance parameter for causality, Hard enforcement, 1e-4,Soft enforcement, 1e-3
OP.EC_PULSE_TOL = xls_parameter(parameter, 'Enforce Causality pulse start tolerance', false, 0.01); % Tolerance parameter for causality, Hard enforcement, 0.05, Soft enforcement, .01
OP.pkg_len_select = xls_parameter(parameter, 'z_p select', true, 1);  % List of package length indexes used to run COM
OP.RX_CALIBRATION = xls_parameter(parameter, 'RX_CALIBRATION', false, false); % Turn on RX_Calibration loop
OP.sigma_bn_STEP = xls_parameter(parameter, 'Sigma BBN step', false, 5e-3); % BBN step for Rx Calibration in volts. Defaults is 0.5e-3
OP.BBN_Q_factor = xls_parameter(parameter, 'BBN Q factor', false, 5);   % Overrides NEXT/FEXT noise Qfactor for  'Force BBN Q factor' used for reporting. does not affect COM.
OP.force_BBN_Q_factor = xls_parameter(parameter, 'Force BBN Q factor', false, false); % Used for reporting and bathtub curves. does not affect COM.
OP.transmitter_transition_time = xls_parameter(parameter, 'T_r', true , 8e-3); % 20% to 80% transition time used for the Gaussian shaped source
OP.RL_norm_test=xls_parameter(parameter, 'ERL_FOM', false, 1); % Defaults to 1 indicating variance is used for FOM determination.  Do not change.

OP.T_r_meas_point = xls_parameter(parameter, 'T_r_meas_point', false, 0); % included for earlier version support. Not recommended to use.
OP.T_r_filter_type= xls_parameter(parameter, 'T_r_filter_type', false, 0);% included for earlier version support. Not recommended to use.
OP.FORCE_TR = xls_parameter(parameter, 'FORCE_TR', false, false);% Included for earlier version support but should be set to 1 in most later config sheets.
% Control with OP.T_r_filter_type and OP.T_r_meas_point for backward
% compatibility
if OP.FORCE_TR
    OP.T_r_meas_point=0;
    OP.T_r_filter_type=1;
end
OP.TDR = xls_parameter(parameter, 'TDR', false, false); % Set to 1 to produce TDR results
OP.TDR_duration= xls_parameter(parameter, 'TDR_duration', false, 5); % only used if N*UI is longer than the TDR duration time.  Default is 5 times the raw s-parameter transit time.
OP.N = xls_parameter(parameter, 'N', false, 0); % duration time in UI which is used for ERL (PTDR)
OP.WC_PORTZ = xls_parameter(parameter, 'WC_PORTZ', false, false); % Do not use: Obsolete. 
OP.T_k= xls_parameter(parameter, 'T_k', false, .6)*1e-9; % Time span (ns) for which the impedance of port is determined using TDR.
OP.ERL_ONLY = xls_parameter(parameter, 'ERL_ONLY', false,0); % Compute ERL only
OP.ERL=xls_parameter(parameter, 'ERL', false, false); % Enables ERL. Needs TDR to be set as well.
if OP.ERL
    OP.PTDR=1;
else
    OP.PTDR=0;
end % ERL needs to do a TDR
OP.SHOW_BRD= xls_parameter(parameter, 'SHOW_BRD', false,0);% indclude added board (PCB) in TDR and ERL. Default is 0.
if OP.WC_PORTZ , OP.TDR=1;end % Obsolete: WC_PORTZ needs to do a TDR
OP.TDR_W_TXPKG = xls_parameter(parameter, 'TDR_W_TXPKG', false,0);% adds tx package for TDR, PTDR, and ERL. Default is 0.
OP.Bessel_Thomson=xls_parameter(parameter, 'Bessel_Thomson', false, false); % enable Bessel Thomsen filter for COM
OP.TDR_Butterworth=xls_parameter(parameter, 'TDR_Butterworth', false, true); % enable Butterworth filter for TDR, PTDR, and ERL
OP.Butterworth=xls_parameter(parameter, 'Butterworth', false, 1); % Enable Butterworth Rx filter for COM compuatetopm
OP.Raised_Cosine=xls_parameter(parameter, 'Raised_Cosine', false,0); % Not used if 0. Default is zero. Should set BT and BW to false
OP.inc_reflect_board=xls_parameter(parameter, 'inc_reflect_board', false,0); % Not used if 0. Default is zero.
OP.AUTO_TFX=xls_parameter(parameter, 'AUTO_TFX', false,0); % Mostly used for device ERL. If sent to 1 the fixture tfx will be estimated.
OP.LIMIT_JITTER_CONTRIB_TO_DFE_SPAN = xls_parameter(parameter, 'LIMIT_JITTER_CONTRIB_TO_DFE_SPAN', false, false); % Experimental. Default is 0.
%OP.impulse_response_truncation_threshold = xls_parameter(parameter, 'Impulse response truncatio threshold', false, 1e-3);
OP.impulse_response_truncation_threshold = xls_parameter(parameter, 'Impulse response truncation threshold', false, 1e-3); % zero padding threshold in fraction of IR peak for the impulse response. Effectively controls the length of time for the PR. Larger values decrease run time and accuracy. Default is 1e-3.
OP.interp_sparam_mag = xls_parameter(parameter, 'S-parameter magnitude extrapolation policy', false, 'linear_trend_to_DC'); % magnitued extrapolation method
OP.interp_sparam_phase = xls_parameter(parameter, 'S-parameter phase extrapolation policy', false, 'extrap_cubic_to_dc_linear_to_inf'); % phase extrapolation method
OP.PMD_type= xls_parameter(parameter, 'PMD_type', false,'C2C'); %  Either C2C or C2M. C2M is for computing VEC and VEO
OP.PHY= xls_parameter(parameter, 'PHY', false, OP.PMD_type); % The keyword OP.PMD_type is now used
if strcmpi(OP.PHY,'C2M') 
    OP.EW=true;
else
    param.T_O=0; % make sure when c2c that sample is at Ts
end
if param.Min_VEO ~=0 && strcmpi(OP.PHY,'C2C')
    OP.PHY='C2Mcom';
end
OP.TDECQ=xls_parameter(parameter, 'TDECQ', false, 0); % Experimental, for only option is none (0) or vma. Default is 0.
switch lower(OP.TDECQ)
    case {false 'none' 'vma'}
    otherwise
        error('%s unrecognized TDECQ keyword',OP.TDECQ)
end
OP.RUNTAG = xls_parameter(parameter, 'RUNTAG', false, ''); % This string is appended to the begining of results files 
if isnan(OP.RUNTAG), OP.RUNTAG='';end
if isnumeric(OP.RUNTAG), OP.RUNTAG=num2str(OP.RUNTAG);end
OP.CDR=xls_parameter(parameter, 'CDR', false, 'MM');% 12/21 from Yuchun Lu to accomdate 'Mod-MM', Defautt is 'MM'
OP.Optimize_loop_speed_up =xls_parameter(parameter, 'Optimize_loop_speed_up', true , 0);% If set to 0 (or default) normal looping, If set to 1 loop speedup by slightly reducing PD Fbin and FIR_threshold for optimize looping only  
% Parameters for error burst probability calculation. Not officially used
OP.use_simple_EP_model = xls_parameter(parameter, 'Use simple error propagation model', false, false);% Use to calculate burst error rate (not normally used
OP.nburst = xls_parameter(parameter, 'Max burst length calculated', false, 0); % Use to calculate burst error rate (not normally used)
OP.COM_EP_margin = xls_parameter(parameter, 'Error propagation COM margin', false, 0); % Use to calculate  error propogation (not normally used)
OP.USE_ETA0_PSD = xls_parameter(parameter, 'USE_ETA0_PSD', false, 0); % Used eta_0 PSD equaiton for sigma_n. Default is 0. Do not use.
OP.SAVE_CONFIG2MAT = xls_parameter(parameter, 'SAVE_CONFIG2MAT', false, 0); % If set to 1 (default) saves parameters in mat file. Requires DIAGNOSTICS to be set.
OP.PLOT_CM = xls_parameter(parameter, 'PLOT_CM', false, 0); % Display CM plots if set to 1. Default is 0.
OP.fraction_of_F_range_start_extrap_from= xls_parameter(parameter, 'fraction_of_F_range_start_extrap_from', true, 0.75); % Frequency (fb) where high frequency extropolation begins for computing IR. Helps control Gibbs phenomena. defualt is 0.75.
OP.COMPUTE_RILN = xls_parameter(parameter, 'COMPUTE_RILN', false, 0); % Computes RILN default is 0.  FOM_RILN reported
OP.COMPUTE_TDILN = xls_parameter(parameter, 'COMPUTE_TDILN', false, OP.COMPUTE_RILN); %  computes TD ILN from complex freq IL fit. FOM_TDILN reported.
OP.SAVE_KEYWORD_FILE = xls_parameter(parameter, 'SAVE_KEYWORD_FILE', false, 0); % Save csv file of COM parameter (OP) and keywords. Not implemented. 
OP.SNR_TXwC0 = xls_parameter(parameter, 'SNR_TXwC0', false, 0); % Adjust SNR_TX with C0
OP.MLSE = xls_parameter(parameter, 'MLSE', false, 0); % MLSE keyword
OP.RXFFE_FLOAT_CTL =  xls_parameter(parameter, 'RXFFE FLOAT CTL', false, 'Taps'); % select taps (taps) or pulse response (ISI) for floating taps
OP.RXFFE_TAP_CONSTRAINT =xls_parameter(parameter, 'RXFFE TAP CONSTRAINT', false, 'Unity Cursor'); % "Unity sum taps", "Unity Cursor", of unbounded 
if OP.MLSE && param.ndfe==0
        error('At least DFE 1 must be set to use MLSE');
end
OP.TIME_AXIS = xls_parameter(parameter, 'TIME_AXIS', false, 'UI'); % if0 OP.display set pulse response xaxis to seconds or UI
% MNSE parameters
OP.Do_XT_Noise= xls_parameter(parameter, 'Do_XT_Noise', false, 1);
OP.FFE_SNR= xls_parameter(parameter, 'FFE_SNR', false, 1);
OP.Do_Colored_Noise= xls_parameter(parameter, 'Do_Colored_Noise', false, 1);
OP.Do_White_Noise=xls_parameter(parameter, 'Do_White_Noise', false, 0);
OP.FFE_OPT_METHOD=xls_parameter(parameter,'FFE_OPT_METHOD',false,'FV-LMS'); % 'MMSE','FV-LMS', 'WIENER-HOPF', 

% need to make sure TD mode does not invoke FD operations
if OP.TDMODE % need to set GET_FD false of TDMODE
    OP.GET_FD=false;
    OP.ERL_ONLY=0;
    OP.ERL=0;
    OP.PTDR=0;
    OP.TDR=0;
    OP.RX_CALIBRATION=0;
end
if OP.SAVE_CONFIG2MAT || OP.CONFIG2MAT_ONLY
    save(matcongfile ,'parameter');
end


%% At the very end of Parameter reading, swap in the proper Tx and Rx values for package parameters based on pkg name
if ~isempty(param.PKG_NAME)
    if length(param.PKG_NAME) == 1
        param.PKG_NAME = [param.PKG_NAME param.PKG_NAME];
    end
    tx_rx_fields = {'C_pkg_board' 'R_diepad'};
    tx_rx_fields_matrix = {'pkg_Z_c'};
    tx_fields = {'z_p_tx_cases' 'z_p_next_cases' 'z_p_fext_cases' 'pkg_gamma0_a1_a2' 'pkg_tau' 'a_thru' 'a_fext'};
    rx_fields = {'z_p_rx_cases' 'a_next'};
    tx_pkg_name=param.PKG_NAME{1};
    rx_pkg_name=param.PKG_NAME{2};
    tx_pkg_struct=param.PKG.(tx_pkg_name);
    rx_pkg_struct=param.PKG.(rx_pkg_name);
    
    %tx_rx_fields: put the value from the tx package in the Tx position and the value from the rx package in the RX position
    for j=1:length(tx_rx_fields)
        tx_val = tx_pkg_struct.(tx_rx_fields{j});
        rx_val = rx_pkg_struct.(tx_rx_fields{j});
        param.(tx_rx_fields{j}) = [tx_val(1) rx_val(2)];
    end
    
    %tx_rx_fields_matrix:  same as tx_rx_fields but in matrix form
    for j=1:length(tx_rx_fields_matrix)
        tx_val = tx_pkg_struct.(tx_rx_fields_matrix{j})(1,:);
        rx_val = rx_pkg_struct.(tx_rx_fields_matrix{j})(2,:);
        param.(tx_rx_fields_matrix{j}) = [tx_val; rx_val];
    end
    
    %tx_fields:  use only the tx package values
    for j=1:length(tx_fields)
        param.(tx_fields{j}) = tx_pkg_struct.(tx_fields{j});
    end
    
    %rx_fields:  use only the rx package values
    for j=1:length(rx_fields)
        param.(rx_fields{j}) = rx_pkg_struct.(rx_fields{j});
    end

end


%%
function [data, SDD, SDC] = read_p2_s2params(infile, plot_ini_s_params, plot_dif_s_params, ports,OP)
%% FUNCTION :: read_sp4_sparams
%
% Description
%   Read the fid of single-ended 4-port complex S-parameters
%   in Touchstone format 'file' and convert to the internal
%   format using the port transform 'ports'
%
%   Created by Mike Y. He
%   April 22, 2005
%
%   Reused some code from
%   Anthony Sanders, Alex Deas, Bob Davidov (24 January 2005)
%   for touchstone 4-port S-matrix import.
%
%   Modified (2012-July-27) by Ken Young to match current indexing scheme and
%   optimized for quicker parameter matching and parsing. also, separated out
%   the plotting algorithms into their own sub-function routines
%
%   Modified December 2021 to use read_Nport_touchstone
%   This is faster reader that is capable of reading touchstone with any number of ports
%
% Input Variables (required)
%   infile              -- The s4p file to be read and converted
%   plot_ini_sparams    -- Plot the initial s-parameter information. For debugging purposes
%   plot_dif_s_params   -- Plot the differential s-parameter information. For debugging purposes
%   ports               -- Re-order the port layout
%
% Output/Return Variables
%   data        -- structure containing network parameter data points and frequency axis
%   sdc         -- the differential in/common-mode out s-parameter data matrix
%   sdd         -- the differential in/differential out s-parameter data matrix
%


% backwards compatibility settings. can be removed in updated code.
if ~exist('OP', 'var'); OP.DISPLAY_WINDOW = true; end
if isempty(ports); ports = [1 2]; end % default order normally used.
ports = [1 2];


if OP.DISPLAY_WINDOW
    set(0,'defaulttextinterpreter','none') % prevents subscripting character in displayed messages
    hMsgBox = msgbox(infile, 'Reading S-Parameter File'); % display a progress bar for reading the s-parameter file(s)
end

%AJG:  fast touchstone read for any number of ports
[sch,schFreqAxis]=read_Nport_touchstone(infile,ports);



D=NaN(size(sch));
% calculate differential s parameter matrix from single ended
for i=1:size(sch,1)
    S(:,:) = sch(i,:,:);
    T = [1 1  ; 1 -1  ];
    W = T * (S / T);
    D(i,:,:) = W(:,:);
end

% D matrix should be
% Scc11 Scd11 Scc12 Scd21
% Sdc11 Sdd11 Sdc12 Sdd12
% Scc21 Scd21 Scc22 Scd22
% Sdc21 Sdd21 Sdc22 Sdd22

% proper values
%AJG:  matrix can be properly referenced after fixing mapping
SDD(:,1,1) = D(:,2,2);
SDC(:,1,1)=  D(:,2,1);
SCC(:,1,1)=  D(:,1,1);
SCD(:,1,1)=  D(:,1,2);



% backwards compatibility output variables
data.m = sch;
%schFreqAxis=schFreqAxis(1:freqCounter); % truncating preallocated array to number of freq points.
data.freq = schFreqAxis;
colors = 'rgbk';

if (plot_ini_s_params == 1)
    figure('name', 'Single-ended s-parameters');set(gcf,'Tag','COM');
    for mj=1:4
        %         subplot(2,2,mj);
        for mi=1:4
            plot(data.freq, 20*log10(abs(data.m(:,mj,mi)+1.0e-15)), ...
                colors(mi), 'linewidth', 2, 'disp', sprintf('S%d%d', mj, mi));
            hold on
        end
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        legend show
        grid on
        title(sprintf('Output port %d', mj));
    end
end
plot_dif_s_params =0;
if (plot_dif_s_params == 1)
    figure('name', 'Mixed-mode s-parameters');set(gcf,'Tag','COM');
    %     subplot(2,1,1);
    for mj=1:1
        for mi=1:1
            plot(data.freq, 20*log10(abs(squeeze(SDD(:,mj,mi)))), ...
                colors((mj-1)*2+mi), 'linewidth',2, 'disp', sprintf('SDD%d%d', mj, mi));
            hold on
        end
    end
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend show
    grid on
    title(infile);
    
    %     subplot(2,1,2);
    %     for mj=1:2
    %         for mi=1:2
    %             plot(data.freq, 20*log10(abs(SDC(:,mj,mi))+1.0e-15), ...
    %                 colors((mj-1)*2+mi), 'linewidth',2, 'disp', sprintf('SDC%d%d', mj, mi));
    %             hold on
    %         end
    %     end
    %     xlabel('Frequency (Hz)');
    %     ylabel('Magnitude (dB)');
    %     legend show
    %     grid on
end

if OP.DISPLAY_WINDOW, close(hMsgBox); end
% end read_sp2_sparam

function [data, SDD, SDC, SCC ] = read_p4_s4params(infile, plot_ini_s_params, plot_dif_s_params, ports,OP,param)
%% FUNCTION :: read_sp4_sparams
%
% Description
%   Read the fid of single-ended 4-port complex S-parameters
%   in Touchstone format 'file' and convert to the internal
%   format using the port transform 'ports'
%
%   Created by Mike Y. He
%   April 22, 2005
%
%   Reused some code from
%   Anthony Sanders, Alex Deas, Bob Davidov (24 January 2005)
%   for touchstone 4-port S-matrix import.
%
%   Modified (2012-July-27) by Ken Young to match current indexing scheme and
%   optimized for quicker parameter matching and parsing. also, separated out
%   the plotting algorithms into their own sub-function routines
%
%  Modified December 2021 to use read_Nport_touchstone
%  This is faster reader that is capable of reading touchstone with any number of ports
%
% Input Variables (required)
%   infile              -- The s4p file to be read and converted
%   plot_ini_sparams    -- Plot the initial s-parameter information. For debugging purposes
%   plot_dif_s_params   -- Plot the differential s-parameter information. For debugging purposes
%   ports               -- Re-order the port layout
%   OP
%   param
% Output/Return Variables
%   data        -- structure containing network parameter data points and frequency axis
%   sdd         -- the differential in/differential out s-parameter data matrix
%   sdc         -- the differential in/common-mode out s-parameter data matrix
%   scc         -- the common mode in/common-mode out s-parameter data matrix
%
%


% backwards compatibility settings. can be removed in updated code.
if ~exist('OP', 'var'); OP.DISPLAY_WINDOW = true; end
if isempty(ports); ports = [1 3 2 4]; end % default order normally used.

% adjust ports to maintain the meaning [in1, in2 , out1, out2] when one
% pair is reversed.
%ports_adj=ports; for k=1:4, ports_adj(k)=find(ports==k); end; ports=ports_adj;

if OP.DISPLAY_WINDOW
    hMsgBox = msgbox(infile, 'Reading S-Parameter File'); % display a progress bar for reading the s-parameter file(s)
end

%AJG:  fast touchstone read for any number of ports
[sch,schFreqAxis]=read_Nport_touchstone(infile,ports);
% matrix to introduce p or n skew on Tx or Rx RIM 12/29/2023
% Sigma's will be form exp(2i*pi*f*skew*1e-12). i.e. if skew = 0 sigma = 1
% need to swap sigma for 1 and 3 and 2 and 4 not RIM 12/29/2023
Sigfct = ...
    @(sigma2,sigma1,sigma4,sigma3)reshape([sigma1.^2,sigma1.*sigma2,sigma1.*sigma3,sigma1.*sigma4,sigma1.*sigma2,sigma2.^2,sigma2.*sigma3,sigma2.*sigma4,sigma1.*sigma3,sigma2.*sigma3,sigma3.^2,sigma3.*sigma4,sigma1.*sigma4,sigma2.*sigma4,sigma3.*sigma4,sigma4.^2],[4,4]);
D=NaN(size(sch));
% calculate differential s parameter matrix from single ended
% skew added RIM 12/29/2023
for i=1:size(sch,1)
    f=schFreqAxis(i);
    sigma_matrix=Sigfct(exp(2i*pi*f*param.Txpskew*1e-12),exp(2i*pi*f*param.Txnskew*1e-12),exp(2i*pi*f*param.Rxpskew*1e-12),exp(2i*pi*f*param.Rxnskew*1e-12) );
    S(:,:) = sch(i,:,:);
    Snew=sigma_matrix.*S;
    T = [1 1 0 0 ; 1 -1 0 0 ; 0 0 1 1 ; 0 0 1 -1];
    W = T * (Snew / T);
    D(i,:,:) = W(:,:);
end

% D matrix should be
% Scc11 Scd11 Scc12 Scd21
% Sdc11 Sdd11 Sdc12 Sdd12
% Scc21 Scd21 Scc22 Scd22
% Sdc21 Sdd21 Sdc22 Sdd22

% proper values
SDD(:,1,1) = D(:,2,2);
SDD(:,2,2) = D(:,4,4);
SDD(:,1,2) = D(:,2,4);
SDD(:,2,1) = D(:,4,2);

SDC(:,1,1) = D(:,2,1);
SDC(:,2,2) = D(:,4,3);
SDC(:,1,2) = D(:,2,3);
SDC(:,2,1) = D(:,4,1);

SCC(:,1,1) = D(:,1,1);
SCC(:,2,2) = D(:,3,3);
SCC(:,1,2) = D(:,1,3);
SCC(:,2,1) = D(:,3,1);

% backwards compatibility output variables
data.m = sch;
%schFreqAxis=schFreqAxis(1:freqCounter); % truncating preallocated array to number of freq points.
data.freq = schFreqAxis;
colors = 'rgbk';

if (plot_ini_s_params == 1)
    figure('name', 'Single-ended s-parameters');set(gcf,'Tag','COM');
    for mj=1:4
        subplot(2,2,mj);
        for mi=1:4
            plot(data.freq, 20*log10(abs(data.m(:,mj,mi)+1.0e-15)), ...
                colors(mi), 'linewidth', 2, 'disp', sprintf('S%d%d', mj, mi));
            hold on
        end
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        legend show
        grid on
        title(sprintf('Output port %d', mj));
    end
end
plot_dif_s_params =0;
if (plot_dif_s_params == 1)
    figure('name', 'Mixed-mode s-parameters');set(gcf,'Tag','COM');
    %     subplot(2,1,1);
    for mj=1:2
        for mi=1:2
            plot(data.freq, 20*log10(abs(SDD(:,mj,mi))), ...
                colors((mj-1)*2+mi), 'linewidth',2, 'disp', sprintf('SDD%d%d', mj, mi));
            hold on
        end
    end
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend show
    grid on
    title(infile);
    %
    %     subplot(2,1,2);
    %     for mj=1:2
    %         for mi=1:2
    %             plot(data.freq, 20*log10(abs(SDC(:,mj,mi))+1.0e-15), ...
    %                 colors((mj-1)*2+mi), 'linewidth',2, 'disp', sprintf('SDC%d%d', mj, mi));
    %             hold on
    %         end
    %     end
    %     xlabel('Frequency (Hz)');
    %     ylabel('Magnitude (dB)');
    %     legend show
    %     grid on
end

if OP.DISPLAY_WINDOW, close(hMsgBox); end
% end read_sp4_sparam
function param_struct = read_package_parameters(parameter,param_struct)

%With the introduction of Package sections in COM spreadsheet, it makes sense to have a single function that grabs all package parameters from a parameter block
%This block should eventually replace what is in read_ParamConfigFile
%It can be called as:  param = read_package_parameters(parameter, param)

if nargin<2
    %param_struct doesn't need to be passed when building a new package structure
    %it is only needed when appending to regular param structure
    param_struct=struct;
end

param_struct.C_pkg_board = xls_parameter(parameter, 'C_p', true)*1e-9; % C_p in nF (single sided)
param_struct.R_diepad = xls_parameter(parameter, 'R_d', true); % Die source termination resistance  (single sided)

param_struct.a_thru = xls_parameter(parameter, 'A_v', true); % Victim differential peak source output voltage (half of peak to peak)
param_struct.a_fext = xls_parameter(parameter, 'A_fe', true); % FEXT aggressor differential peak source output voltage (half of peak to peak)
param_struct.a_next = xls_parameter(parameter, 'A_ne', true); % NEXT aggressor differential peak source output voltage (half of peak to peak)

param_struct.z_p_tx_cases = xls_parameter(parameter, 'z_p (TX)', true).'; % List of victim transmitter package trace lengths in mm, one per case
[ncases, mele]=size(param_struct.z_p_tx_cases);
if mele ==2
    param_struct.flex=2;
elseif mele==4
    param_struct.flex=4;
elseif mele==1
    param_struct.flex=1;
else
    error('config file syntax error')
end
param_struct.z_p_next_cases = xls_parameter(parameter, 'z_p (NEXT)', true).'; % List of NEXT transmitter package trace lengths in mm, one per case
[ncases1, mele1]=size(param_struct.z_p_next_cases);
if ncases ~= ncases1 || mele ~= mele1
    error('All TX, NEXT, FEXT, Rx cases must agree');
else
end
param_struct.z_p_fext_cases = xls_parameter(parameter, 'z_p (FEXT)', true).'; % List of FEXT transmitter package trace lengths in mm, one per case
[ncases1, mele1]=size(param_struct.z_p_fext_cases);
if ncases ~= ncases1 ||  mele ~= mele1
    error('All TX, NEXT, FEXT, Rx cases must agree');
else
end
param_struct.z_p_rx_cases = xls_parameter(parameter, 'z_p (RX)', true).'; % List of FEXT receiver package trace lengths in mm, one per case
[ncases1, mele1]=size(param_struct.z_p_rx_cases);
if ncases ~= ncases1 ||  mele ~= mele1
    error('All TX, NEXT, FEXT, Rx cases must agree');
else
end
% Table 93A-3 parameters
param_struct.pkg_gamma0_a1_a2 = xls_parameter(parameter, 'package_tl_gamma0_a1_a2', true, [0 1.734e-3 1.455e-4]); %Fitting parameters for package model per unit length. First element is in 1/mm and affects DC loss of package model . Second element is in ns1/2/mm and affects loss proportional to sqrt(f). Third element is in ns/mm and affects loss proportional to f.
param_struct.pkg_tau = xls_parameter(parameter, 'package_tl_tau', true, 6.141e-3); % Package model transmission line delay ns/mm
param_struct.pkg_Z_c = xls_parameter(parameter, 'package_Z_c', true, 78.2).';% Package model transmission line characteristic impedance [ Tx , Rx ]
[ ncases1, mele1]=size(param_struct.pkg_Z_c);%
if   mele ~= mele1
    error('tx rx pairs must have thesame number element entries as TX, NEXT, FEXT, Rx');
else
end
if mele1==2 % fuill in a array if only a 2 element flex package is specified
    for ii=1:ncases
        param_struct.z_p_fext_casesx(ii,:)=  [param_struct.z_p_fext_cases(ii,:)' ;[ 0 ; 0 ]]';
        param_struct.z_p_next_casesx(ii,:)=  [param_struct.z_p_next_cases(ii,:)' ;[ 0 ; 0 ]]';
        param_struct.z_p_tx_casesx(ii,:)=  [param_struct.z_p_tx_cases(ii,:)' ;[ 0 ; 0 ]]';
        param_struct.z_p_rx_casesx(ii,:)=  [param_struct.z_p_rx_cases(ii,:)' ;[ 0 ; 0 ]]';
    end
    param_struct.z_p_fext_cases =  param_struct.z_p_fext_casesx;
    param_struct.z_p_next_cases=  param_struct.z_p_next_casesx;
    param_struct.z_p_tx_cases=  param_struct.z_p_tx_casesx;
    param_struct.z_p_rx_cases=  param_struct.z_p_rx_casesx;
    param_struct.pkg_Z_c=[param_struct.pkg_Z_c' ;[ 100 100 ; 100 100 ]]';
end
function [chdata,SDDch,SDDp2p] = read_s4p_files(param, OP, chdata)
%% extract s-parameter and convert to differential mode
% extract s-parameter data from files and apply tx and rx filters as well as package filters
num_files=length(chdata);
if ~OP.DISPLAY_WINDOW, fprintf('reading file '); end
for i=1:num_files
    if OP.DISPLAY_WINDOW; hwaitbar=waitbar(0);end
    progress = i/num_files;
    if OP.DISPLAY_WINDOW
        [~,a]=fileparts(chdata(i).filename);
        waitbar(progress, hwaitbar, ['Processing ' a]); figure(hwaitbar); drawnow;
    else
        fprintf('%i ',i);
    end
    
    % Skip reading file if it was already read (multiple test cases)
    if (~isfield(chdata(i), 'faxis')) || isempty(chdata(i).faxis)
        switch lower(chdata(i).ext)
            case '.s2p' % for differential return loss
                [Sch,SDDch] = read_p2_s2params(chdata(i).filename,  0, 0, param.snpPortsOrder, OP);
                chdata(i).fmaxi = length(Sch.freq);
                chdata(i).faxis = Sch.freq;
                chdata(i).sdd11_raw = transpose(SDDch(1:chdata(i).fmaxi,1,1));
                SDDp2p(i)=NaN;
                chdata(i).sdd11_orig=chdata(i).sdd11_raw;
                chdata(i).sdd11=chdata(i).sdd11_raw;
            case '.s4p'
                if length(param.snpPortsOrder) ~= 4
                    error( 'warning:sNpFilePortMismatch', ...
                        '\n\t The number of ports defined (%G) does not match the sNp file type (%s)', ...
                        length(param.snpPortsOrder), ...
                        chdata(i).ext ...
                        );
                end
                % read function returns differnetial mode parameters
                if param.package_testcase_i==1 % added to speed up cases e.g. don't read file in twice
                    [Sch, SDDch, SDCch] = read_p4_s4params(chdata(i).filename,  0, 0, param.snpPortsOrder, OP,param);
%                     param.holdsdata(i).Sch= Sch;
%                     param.holdsdata(i).SDDch=  SDDch;
%                     param.holdsdata(i).SDCch= SDCch;
                else
                    error('If this line is reached, there is a logic error');
%                     Sch=param.holdsdata(i).Sch;
%                     SDDch=param.holdsdata(i).SDDch;
%                     SDCch=param.holdsdata(i).SDCch;
                end
                chdata(i).fmaxi = length(Sch.freq);
                

                if Sch.freq(chdata(i).fmaxi) < param.fb
                    warning('COM:read_s4p:MaxFreqTooLow', ...
                        'In %s: the maximum frequency provided, %g, is less than the signaling rate: %g', ...
                        chdata(i).filename, Sch.freq(end), param.fb);
                end
                if Sch.freq(1) > param.max_start_freq
                    warning('COM:read_s4p:StartFreqTooHigh', ...
                        'In %s: minimum frequency, %.2g GHz, is larger than the recommended %.2g GHz', ...
                        chdata(i).filename, Sch.freq(1)/1e9, param.max_start_freq/1e9);
                end
                freqstep=diff(Sch.freq);
                % ignore frequency differences up to 1 Hz - possible numerical artifacts
                if max(freqstep)-min(freqstep) > 1
                    warning('COM:read_s4p:NonUniformFreqSpacing', 'In %s: non-uniform frequency steps: min=%.3g GHz, max=%.3g GHz', ...
                        chdata(i).filename, min(freqstep)/1e9, max(freqstep)/1e9);
                end
                if max(freqstep) - param.max_freq_step > 1
                    warning('COM:read_s4p:FreqStepTooHigh', 'In %s: frequency step, %.2g GHz, is larger than the recommended %.2g GHz', ...
                        chdata(i).filename, max(freqstep)/1e9, param.max_freq_step/1e9);
                end
                
                chdata(i).faxis = Sch.freq;
                chdata(i).sdd12_raw = transpose(SDDch(1:chdata(i).fmaxi,1,2));
                chdata(i).sdd21_raw = transpose(SDDch(1:chdata(i).fmaxi,2,1));
                chdata(i).sdd22_raw = transpose(SDDch(1:chdata(i).fmaxi,2,2));
                chdata(i).sdd11_raw = transpose(SDDch(1:chdata(i).fmaxi,1,1));
                % mode conversion
                chdata(i).sdc12_raw = transpose(SDCch(1:chdata(i).fmaxi,1,2));
                chdata(i).sdc21_raw = transpose(SDCch(1:chdata(i).fmaxi,2,1));
                chdata(i).sdc22_raw = transpose(SDCch(1:chdata(i).fmaxi,2,2));
                chdata(i).sdc11_raw = transpose(SDCch(1:chdata(i).fmaxi,1,1));
                %save original and add board (if required)
                chdata(i).sdd11_orig=chdata(i).sdd11_raw;
                chdata(i).sdd22_orig=chdata(i).sdd22_raw;
                chdata(i).sdd12_orig=chdata(i).sdd12_raw;
                chdata(i).sdd21_orig=chdata(i).sdd21_raw;
                if OP.include_pcb
                    % add boards to sdd
                    [chdata(i).sdd11_raw, chdata(i).sdd12_raw, chdata(i).sdd21_raw, chdata(i).sdd22_raw] = add_brd(chdata(i), param, OP);
                    
                end
                %save final return loss (after the boards were included)
                chdata(i).sdd11=chdata(i).sdd11_raw;
                chdata(i).sdd22=chdata(i).sdd22_raw;
            otherwise
                error('Extension "%s" in file "%s" is not supported',chdata(i).ext,chdata(i).filename);
        end
        
        %Crosstalk frequency axis must be the same as Thru
        if i>1
            %error on length difference
            if length(chdata(i).faxis)~=length(chdata(1).faxis)
                error('Crosstalk file "%s" has different number of frequency points',chdata(i).filename);
            end
            %error if any value > 1Hz (don't want to check for exact
            %equality in case of floating point error)
            Fdiff=abs(chdata(i).faxis-chdata(1).faxis);
            if max(Fdiff)>1
                error('Crosstalk file "%s" has a different frequency axis',chdata(i).filename);
            end
        end
    else
        SDDch(:,1,2)=chdata(i).sdd12_raw;
        SDDch(:,2,1)=chdata(i).sdd21_raw;
        SDDch(:,1,1)=chdata(i).sdd11_raw;
        SDDch(:,2,2)=chdata(i).sdd22_raw;
    end
    chdata(i).sigma_ACCM_at_tp0=0;
    if ~param.FLAG.S2P
        if  OP.INC_PACKAGE ~= 0 || (OP.RX_CALIBRATION == 1  && i==1)
            if (OP.RX_CALIBRATION == 1  && i==2)
                chdata(i).sdd21=chdata(i).sdd21_raw;
            else
                %updated package construction with single function for both DD and DC
                [chdata(i).sdd21p,SDDp2p(i)]= s21_pkg(chdata(i), param, OP, i);
                [chdata(i).sdd21p_nodie]= s21_pkg(chdata(i), param, OP, i, 'dd', 0);
                chdata(i).sdd21=chdata(i).sdd21p;
                if 1 % for AC CM noise inclusion
                    [chdata(i).sdc21p,SDCp2p(i),chdata(i).sigma_ACCM_at_tp0]= s21_pkg(chdata(i), param, OP, i,'dc');
                    chdata(i).sdc21=chdata(i).sdc21p;
                end
            end
        else
            chdata(i).sdd21=chdata(i).sdd21_raw;
        end
        chdata(i).sdd21f=chdata(i).sdd21_orig; % used for FD analysis i.e. not filtered (RIM 9/24/2021 without boards or packages)
    end
end
if ~OP.DISPLAY_WINDOW, fprintf('\n'); end

function result = readdataSnPx(filename, nport)
%function [freq, cs] = readdataSnPx(filename, nport)
% [freq, cs] = readdataSnP(filename, nport, format, nheader)
%
% Read Touchstone file with frequencies in units of Hertz
%
% Input:
% ======
% filename: Name of the Touchstone/SnP file
% nport: Number of ports
% format: 'RI' for real/imag, 'MA' for mag/angle (check option line in the
%         Touchstone file)
% nheader: Number of header lines (comment lines plus option line in the
%          Touchstone file)
%
% Output:
% =======
% freq: Vector of frequencies [Hz]
% cs: 3-D array of complex-valued S parameters where cs(i,j,k) is S(i,j)
%     at frequency freq(k)
%
% Note: If frequency unit is not Hertz (but GHz, MHz etc.) simply scale
% frequencies appropriately after reading the data.
%
% Ref.: Touchstone(R) File Format Specification, Rev.1.1,
% EIA/IBIS Open Forum, 2002.
%
% Written by Henning Braunisch, September 2004.
% Updated by Steven Krooswyk, April 2006.


fid = fopen(filename, 'r');


% Skip header lines
str = ' ';
n = 0;
while ~strcmp(str(1),'#')
    str = fgetl(fid);
    if isempty(str)
        str=' ' ;
        if n > 1000
            display('error: could not find config line (#)')
            break
        end
    end
    n = n + 1;
end

% parse configuration line
A=sscanf(str,'%1s %2s %1s %2s %1s %2s',[1,inf]);
p = find(A=='S');           %position of 'S'
units = lower(A(2:p-1));    %units before 'S'
format = A(p+1:p+2);        %format after 'S'

% skip any more header lines
%while ~str

nk = 0; % frequency counter
while 1
    
    [temp, count] = fscanf(fid, '%f', 1);
    if count == 0
        temp2 = fscanf(fid, '%s', 1);
        if ~isempty(temp2), fgetl(fid); continue, end;
        break
    end
    nk = nk+1; freq(1,nk) = temp;  %#ok<AGROW>
    for ni = 1:nport
        for nj = 1:nport
            switch lower(format)
                case 'ma'
                    mag = fscanf(fid, '%f', 1); ang = fscanf(fid, '%f', 1);
                    cs(ni,nj,nk) = mag * exp(1i*ang*pi/180); %#ok<AGROW>
                case 'ri'
                    re = fscanf(fid, '%f', 1); im = fscanf(fid, '%f', 1);
                    cs(ni,nj,nk) = complex(re, im); %#ok<AGROW>
                case 'db'
                    db = fscanf(fid, '%f', 1); ang = fscanf(fid, '%f', 1);
                    M = 10^(db/20);
                    %re = M*cos(ang);
                    %im = M*sin(ang);
                    re =  M*cos(ang * pi / 180);
                    im =  M*sin(ang * pi / 180);
                    cs(ni,nj,nk) = complex(re, im); %#ok<AGROW>
                otherwise
                    error('readdataSnP: Unknown data format');
            end
        end
    end
end

fclose(fid);

% If 2-port then swap S_12 and S_21 per Touchstone spec
if nport == 2
    temp = cs(2,1,:);
    cs(2,1,:) = cs(1,2,:);
    cs(1,2,:) = temp;
end

% Update freq units to Hz
switch lower(units)
    case 'hz'
        
    case 'khz'
        freq=freq.*1e3;
    case 'mhz'
        freq=freq.*1e6;
    case 'ghz'
        freq=freq.*1e9;
end

% passivity check
result.freq = freq;
result.cs    = cs;

function recolor_plots(ax)

if ~verLessThan('matlab', '8.4.0')
    return
end
colors='brgcmk';
ch=flipud(get(ax, 'children'));

for k=1:length(ch)
    set(ch(k), 'Color', colors(mod(k-1, length(colors))+1));
    set(ch(k), 'LineWidth', 2*floor((k-1)/length(colors))+1);
end
legend (ax, 'off');
warning('off', 'MATLAB:legend:PlotEmpty');
set(legend (ax, 'show'), 'interp', 'none');

function result = reduce(var1)
% --- Reduce 1x1xn array to 1xn (aka squeeze)
out = zeros(1,length(var1));
out(1,:) = var1(1,1,:);
result=out;

function [s21p,SCH,sigma_ACCM_at_tp0] = s21_pkg(chdata, param, OP, channel_number,mode,include_die)
% concatenates package reflections with s21,s11,and s22 with spec return loss (gammas)
% faxis is the frequency array
% s21, s11, s22 are the corresponding array of differential parameters
% s21p includes the VFT and Tx filter if include_die=1
if nargin<6
    include_die=1;
end
if nargin<5
    mode='dd';
end

s21=chdata.(['s' mode '21_raw']);
s12=chdata.(['s' mode '12_raw']);
s11=chdata.(['s' mode '11_raw']);
s22=chdata.(['s' mode '22_raw']);
faxis=chdata.faxis;
channel_type=chdata.type;

if strcmpi(mode,'dd')
    s11=s11*param.kappa1;
    s22=s22*param.kappa2;
end


Z0=param.Z0;
%sigma_ACCM_at_tp0 is only used when mode=DC
sigma_ACCM_at_tp0=0;

% The following three parameters have possibly different valuesF for TX and
% RX (so can be 2-element vectors).
R_diepad = param.R_diepad;

%Make TX Package
[ s11out, s12out, s21out, s22out]=make_full_pkg('TX',faxis,param,channel_type,mode,include_die);

%Make RX Package
%Important:  RX pkg doesn't change based on mode being dd or dc.  So 'dd' is always passed
[ s11in, s12in, s21in, s22in]=make_full_pkg('RX',faxis,param,channel_type,'dd',include_die);


%     p(1 ,1, :)=s11in;
%     p(2 ,2, :)=s22in;
%     p(1 ,2, :)=s12in;
%     p(2 ,1, :)=s21in;
%
%     S=sparameters(p,faxis);
%     rfwrite(S,'temp.s4p');

if strcmpi(mode,'dc')
    RTX=R_diepad(param.Tx_rd_sel)/2;
    RRX=R_diepad(param.Rx_rd_sel)/2;
    Z0gamma=Z0/2;
else
    RTX=R_diepad(param.Tx_rd_sel);
    RRX=R_diepad(param.Rx_rd_sel);
    Z0gamma=Z0;
end
if OP.IDEAL_TX_TERM || (OP.RX_CALIBRATION == 1 && channel_number == 2) || OP.include_pcb == 2
    gamma_tx=0;
else
    gamma_tx=(RTX-Z0gamma)/(RTX+Z0gamma);% equation 93A-17
end
if OP.IDEAL_RX_TERM
    gamma_rx=0;
else
    gamma_rx=(RRX-Z0gamma)/(RRX+Z0gamma);% equation 93A-17
end

if OP.INC_PACKAGE==0
    s21p= s21;
    warning('do not use INC_PACKAGE = 0. Instead use package parameters)');
else
    if OP.RX_CALIBRATION == 1 && channel_number == 2
        %   for calibration do not include the transmitter package
        [s11out_rx, s12out_rx, s21out_rx, s22out_rx ] = combines4p( s11, s12, s21, s22, s22out, s12out, s21out, s11out ); %#ok<ASGLU> % s22 is ball side of package
        SCH.Frequencies=faxis;
        SCH.Parameters(1,1,:)=s11out_rx;
        SCH.Parameters(2,2,:)=s22out_rx;
        SCH.Parameters(1,2,:)=s12out_rx;
        SCH.Parameters(2,1,:)=s21out_rx;
        SCH.NumPorts=2;
        SCH.Impedance=100;
        %% Equation 93A-18
        if include_die
            s21p= s21out_rx.*(1-gamma_tx).*(1+gamma_rx)./(1.- s11out_rx.*gamma_tx - s22out_rx.*gamma_rx  -s21out_rx.^2.*gamma_tx.*gamma_rx +s11out_rx.*s22out_rx.*gamma_tx.*gamma_rx);
        else
            s21p=s21out_rx; % if no die we do not want a VTF
        end
    else
        %% Equations 93A-4 to 93A-7
        if ~OP.IDEAL_TX_TERM
            [s11, s12, s21, s22] = combines4p( s11out, s12out, s21out, s22out, s11, s12, s21, s22 ); %#ok<ASGLU>
        end
        H_t=ones(1,length(faxis)); % .3bj compatibility
        if OP.IDEAL_TX_TERM ||  OP.T_r_filter_type == 1
            % for RITT testing with good termination as in some instruments
            % and tx filter when required
            if OP.T_r_filter_type==0
                H_t = exp(-(pi*faxis/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
            else
                tr=OP.transmitter_transition_time;
                f9=faxis/1e9;
                if OP.T_r_meas_point == 1
                    k=1.9466+7.12*sqrt(1-6.51e-3/tr);
                    H_t=105./(f9.^4*(k*tr)^4 - f9.^3*(k*tr)^3*10i - 45*f9.^2*(k*tr)^2 + f9*(k*tr)*105i + 105);
                else
                    H_t = exp( -2*(pi*f9*tr/1.6832).^2 ).*exp(-1j*2*pi*f9*tr*3);
                end
                
            end
        end
        if strcmpi(mode,'dc')
%             H_t=ones(1,length(faxis)); % not sure if we need a H_t or not. Is the CM noise correlated to the edge rate?
        end
        if ~OP.IDEAL_RX_TERM
            [s11, s12, s21, s22] = combines4p( s11, s12, s21, s22, s22in, s21in, s12in, s11in ); %#ok<ASGLU> % s22 is ball side of package
        else
            warning('do not use IDEAL_RX_TERM. Instead hard code package and TR parameters')
        end
        %% Equation 93A-18 and part of 93A-1: Ht fix in V290 identified by Bill Kirkland and Ed Frlan ( s21^2 changed to s12*s21 )
        if include_die
            s21p= H_t.*s21.*(1-gamma_tx).*(1+gamma_rx)./(1.- s11.*gamma_tx - s22.*gamma_rx  -s21.*s12.*gamma_tx.*gamma_rx +s11.*s22.*gamma_tx.*gamma_rx);
        else
            s21p=s21; % if no die we do not want a VTF
        end
    end
    
    if strcmpi(mode,'dc')
        % compute AC_CM_RMS at tp0
        OP.TX_BesselThomson=1; %AC CM is measured in scopes with at BT filter
        H_bt=Bessel_Thomson_Filter(param,faxis,OP.TX_BesselThomson);
        if channel_number == 1
            f_int= faxis( faxis<=param.ACCM_MAX_Freq );
            H_cc= s21out.*(1-gamma_tx)./(1.- s11out.*gamma_tx ).*H_bt;
            sigma_ACCM_at_tp0= sqrt(2*param.AC_CM_RMS_TX^2*sum( abs( H_cc(2:length(f_int)) ).^2 .* diff(f_int))/f_int(end)) ;
            %     S=sparameters(p,faxis);
            %     rfwrite(S,'temp.s4p');
        end
    end
    
    SCH.Frequencies=faxis;
    SCH.Parameters(1,1,:)=s11;
    SCH.Parameters(2,2,:)=s22;
    SCH.Parameters(1,2,:)=s12;
    SCH.Parameters(2,1,:)=s21;
    SCH.NumPorts=2;
    if strcmpi(mode,'dc')
        SCH.Impedance=25;
    else
        SCH.Impedance=100;
    end
    
end
function [voltage, t_base, causality_correction_dB, truncation_dB] = ...
    s21_to_impulse_DC(IL, freq_array, time_step, OP)
% Creates a time-domain impulse response from frequency-domain IL data.
% IL does not need to have DC but a corresponding frequency array
% (freq_array) is required.
%
% Causality is imposed using the Alternating Projections Method. See also:
% Quatieri and Oppenheim, "Iterative Techniques for Minimum Phase Signal
% Reconstruction from Phase or Magnitude", IEEE Trans. ASSP-29, December
% 1981 (http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163714)

ILin=IL;
fmax=1/time_step/2;
freq_step=(freq_array(3)-freq_array(2))/1;
fout=0:1/round(fmax/freq_step)*fmax:fmax;
if all(IL==0)
    %response with all zeros is problematic.  set to all eps and avoid interp function
    IL=ones(1,length(fout))*eps;
else
    IL=interp_Sparam(ILin,freq_array,fout, OP.interp_sparam_mag, OP.interp_sparam_phase,OP);
    IL_nan = find(isnan(IL));
    for in=IL_nan
        IL(in)=IL(in-1);
    end
end
IL = IL(:);
% add padding for time steps
% IL_symmetric = [IL(1:end-1);0; flipud(conj(IL(2:end-1)))];
IL_symmetric = [real(IL(1)); IL(2:end-1); real(IL(end)); flipud(conj(IL(2:end-1)))];
impulse_response = real(ifft(IL_symmetric));
L = length(impulse_response);
t_base = (0:L-1)/(freq_step*L);

original_impulse_response=impulse_response;
% Correct non-causal effects frequently caused by extrapolation of IL
% Assumption: peak of impulse_response is in the first half, i.e. not anti-causal
abs_ir=abs(impulse_response);
a = find(abs_ir(1:L/2) > max(abs_ir(1:L/2))*OP.EC_PULSE_TOL);
start_ind = a(1);

err=inf;
while ~all(impulse_response==0)
    impulse_response(1:start_ind)=0;
    impulse_response(floor(L/2):end)=0;
    IL_modified=abs(IL_symmetric).*exp(1j*angle(fft(impulse_response)));
    ir_modified = real(ifft(IL_modified));
    delta = abs(impulse_response-ir_modified);
    
    err_prev = err;
    err=max(delta)/max(impulse_response);
    if err<OP.EC_REL_TOL || abs(err_prev-err)<OP.EC_DIFF_TOL
        break;
    end
    
    impulse_response=ir_modified;
end

causality_correction_dB=20*log10(norm(impulse_response-original_impulse_response)/norm(impulse_response));

if ~OP.ENFORCE_CAUSALITY
    impulse_response = original_impulse_response;
end
% truncate final samples smaller than 1e-3 of the peak
ir_peak = max(abs(impulse_response));
ir_last  = find(abs(impulse_response)>ir_peak*OP.impulse_response_truncation_threshold, 1, 'last');

voltage = impulse_response(1:ir_last);
t_base = t_base(1:ir_last);

truncation_dB=20*log10(norm(impulse_response(ir_last+1:end))/norm(voltage));

function S =s_for_c2(zref,f,cpad)
% S is 2 port s parameters out
S_Parameters(1,1,:) =  -1i*2*pi.*f*cpad*zref./(2+1i*2*pi.*f*cpad*zref);
S_Parameters(2,2,:) =  -1i*2*pi.*f*cpad*zref./(2+1i*2*pi.*f*cpad*zref);
S_Parameters(2,1,:) = 2./(2+1i*2*pi.*f*cpad*zref);
S_Parameters(1,2,:) = 2./(2+1i*2*pi.*f*cpad*zref);
S=sparameters(S_Parameters,f,zref);

function S =s_for_c4(zref,f,cpad)

S2 = s_for_c2(zref,f,cpad);
S4P=s2_to_s4(S2.Parameters);
S=sparameters(S4P,f,zref);
S.Parameters=snp2smp(S.Parameters,zref,[ 1 3 2 4]);




%%
function [ cmd_str ] = save_cmd_line( config_file, chdata, num_fext,num_next, cli_name )
% save commmend string
% for saving from interactive queries


cmd_str=[ cli_name '('  '''' config_file '''' ',' num2str(num_fext) ', ' num2str(num_next) ',' '''' chdata(1).filename ''''];
for i=1:num_next+num_fext
    cmd_str= [cmd_str  ',' '''' chdata(i+1).filename ''''];
end
cmd_str= [ cmd_str ')'];


%%%%% require the RF tool box
%%
function [ h ] = savefigs( param, OP )

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
        saveas(h(ii), fullfile(OP.RESULT_DIR, [figname '.fig']));
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

%%
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

%%
function pdf_out = scalePDF(pdf,scale_factor)
                pdf_out=pdf;
                pdf_out.Min=floor(pdf.Min*scale_factor);
                pdf_out.x=(pdf_out.Min:-pdf_out.Min)*pdf_out.BinSize;
                pdf_out.y=interp1(pdf.x*scale_factor,pdf.y,pdf_out.x);
                pdf_out.y(1)= pdf_out.y(2); % NAN interp work around
                pdf_out.y(end)= pdf_out.y(end-1); % NAN interp work around
                pdf_out.y=pdf_out.y/sum(pdf_out.y);
function t_params = stot(s_params)
% p 67 R. Mavaddat. (1996). Network scattering parameter. Singapore: World Scientific.
% ISBN 978-981-02-2305-2. http://books.google.com/?id=287g2NkRYxUC&lpg=PA65&dq=T-parameters+&pg=PA67.
[s11, s12, s21, s22] = deal(s_params(1,1,:), s_params(1,2,:), s_params(2,1,:), s_params(2,2,:));
delta = (s11.*s22-s12.*s21);
s21(s21==0)=eps;
t_params = [1./s21, -s22./s21; s11./s21, -delta./s21];

function csv_string = str2csv(c)
% convert a cell array of strings to a csv string
cell_tmp = cell(2, length(c));
cell_tmp(1,:)=c;
cell_tmp(2,:) = {','};
cell_tmp{2,end} = '';
csv_string=strcat(cell_tmp{:});

function [s11, s12, s21, s22] = synth_tline(f, Z_c, Z_0, gamma_coeff, tau, d)
f_GHz=f/1e9;
%% Equation 93A-10 %%
gamma_1 = gamma_coeff(2)*(1+1i);
%% Equation 93A-11 %%
gamma_2 = gamma_coeff(3)*(1-2i/pi*log(f_GHz)) + 2i*pi*tau;
%% Equation 93A-9 %%
gamma = gamma_coeff(1)+gamma_1.*sqrt(f_GHz)+gamma_2.*f_GHz;
gamma(f_GHz==0) = gamma_coeff(1);

%% Equation 93A-12 %%
if d==0
    %force matched impedance if length is 0
    %otherwise divide by zero can occur if Z_c=0
    rho_rl=0;
else
    rho_rl=(Z_c-2*Z_0)/(Z_c+2*Z_0);
end

exp_gamma_d = exp(-d*gamma);
%% Equations 93A-13 and 93A-14 %%
s11 = rho_rl*(1-exp_gamma_d.^2)./(1-rho_rl^2*exp_gamma_d.^2);
s21 = (1-rho_rl^2)*exp_gamma_d./(1-rho_rl^2*exp_gamma_d.^2);
s12 = s21;
s22 = s11;

function s_params = ttos(t_params)
% p 67 R. Mavaddat. (1996). Network scattering parameter. Singapore: World Scientific.
% ISBN 978-981-02-2305-2. http://books.google.com/?id=287g2NkRYxUC&lpg=PA65&dq=T-parameters+&pg=PA67.
[t11, t12, t21, t22] = deal(t_params(1,1,:), t_params(1,2,:), t_params(2,1,:), t_params(2,2,:));
delta = t11.*t22-t21.*t12;
t11(t11==0)=eps;
s_params = [t21./t11, delta./t11; 1./t11, -t12./t11];

function [out_var,varg_out]=varargin_extractor(varargin)

if isempty(varargin)
    out_var=[];
    varg_out={};
else
    out_var=varargin{1};
    varg_out=varargin;
    varg_out(1)=[];
end
    
    
function results= vma(PR, M)
% PR=sbr.Data;
% M=32;
% PR is the pulse response
% M is samples per UI
[ seq, syms, syms_nrz ] = PRBS13Q( );
% seq uses [ -1 -1/3 1/3 1] & syms uses [ 0 1 2 3 ]
symbols=seq; 
imaxPR=find(PR==max(PR),1,'first'); % find index for peak
% start end symbols index for 7 3's and 6 0's
indx_S3x7_start=M*(strfind(syms,ones(1,7)*3)-1)+imaxPR;
indx_S3x7_end=M*(strfind(syms,ones(1,7)*3)+5)+imaxPR;
indx_S0x6_start=M*(strfind(syms,ones(1,6)*0))-1+imaxPR;
indx_S0x6_end=M*(strfind(syms,ones(1,6)*0)+4)+imaxPR;
% superposition code
shifting_vector=kron(symbols,[ 1  zeros(1,M-1) ]) ;
Bit_stream_response=filter(PR,1, shifting_vector);
% find center of 3's and 0's
icent3=floor((indx_S3x7_end-indx_S3x7_start)/2 + indx_S3x7_start);
icent0=floor((indx_S0x6_end-indx_S0x6_start)/2 + indx_S0x6_start);
% plot(Bit_stream_response(indx_S3x7_start:indx_S3x7_end))
% hold on
% plot(Bit_stream_response(indx_S0x6_start:indx_S0x6_end))
P_3 = mean( Bit_stream_response((icent3-M):(icent3+M) ) );
P_0 = mean( Bit_stream_response((icent0-M):(icent0+M) ) );
VMA= P_3 - P_0;
results.P_3=P_3;
results.P_0=P_0;
results.VMA=VMA;
function line_intersection=vref_intersect(eye_contour,x_in,vref)

%slope of the 2 sample points around vref crossing
m1=(eye_contour(x_in,1)-eye_contour(x_in-1,1));
%x-intercept for the line
b1=eye_contour(x_in,1)-m1*x_in;
% drawing a horizontal line through vref so slope = 0
m2=0;
%special case for horizontal line, b=y
b2=vref;
%the x-value of line intersection = (b2-b1)/(m1-m2)
%sinc m2 is always 0 and b2 is always vref, this could be stated as (vref-b1)/m1
%And usually vref is 0, so it further reduces to -b1/m1
line_intersection=(b2-b1)/(m1-m2);





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
function zzz_list_of_changes
% structures:
% chdata(i)
%           i= 1 --> THRU index
%           i= 2, num_fext+1 --> FEXT channel index
%           i= num_fext+2, num_next+num_fext+1
%       base: name of THRU file
%          A: amplitude
%       type: 'THRU', 'NEXT', or 'FEXT'
%        ftr: Rise time frequency
%      fmaxi: max number of frequency points
%      faxis: frequency array [Hz]
%      sdd21: (Htot) rewritten as product of ,vtf based on pkg RL ,tx filter, Rx filter
%      sdd22: differential RL
%      sdd11: differential RL
%     sdd21p: vtf based on pkg RL , this an interim parameter and set to sdd21
%     sdd21f: raw differential IL not filtered use for FD plots
% added output_args.peak_uneq_pulse_mV
% added output_args.cable_loss when "Include PCB" is not 0 in the config file
% added: tap c(-2) c(2) and c(3)
% added: g_DC_HP and f_HP_PZ
% added: new value for "Include PCB" = 2 for cable Rx compliance test only Rx host added
% added: BREAD_CRUMBS is 1 then a mat file with the structures params and OP is created in the results directory
% added: T_r_filter_type for RITT testing when IDEAL_TX_TERM is 1:
% added T_r_meas_point for RITT, if 0, measurement was at tp0, if 1 measurement was tp0a
% 0 is for is for Gaussian filter and 1 is for a 4th order Bessel-Thomson filter
% fixed INCLUDE_CTLE=0 to really remove from computation
% r161a fixed matlab version issue when for OP.INCLUDE_CTLE=0
% r162 adjusting RITT rise time to Mike Dudek's recommendations also always enable risetime filter if T_r_filter_type=1
% r162 tx and rx package impedance {Zc)
% r162a Gaussian equation corrected
% r163 cast snr_tx with package test case
% r164 fix pdf for very low noise and lo pass filter enhancements
% r164 add zero gain at nqyist CTLE as in CL12e
% r165 add simpler congfig command called FORCE_TR (force risetime)
% r200 cm3 and cm4 added cp3 removed
% r200 fixed problem in s21_to_impulse_DC when s parameter have a DC entry
% r200 ILD_FOM updated to EQ93A-55  ERL adde
% r200 improved phase interpolation for return loss time conversion
% r200a db = @(x) 20*log10(abs(x)) added so sig processing toolbox won't be required
% r200b Fixed error in bifurcation of Tx/Rx Rd%
% r200c missed on fix for interpolation
% r210 new ERL with time gating function
% r224 update ERL with from D3.1
% r226  fix s2p reading problem
%    change SNR_ISI_XTK_normalized_1_sigma to SNR_ISI (with nulled noise)
%    Fix Rx calibration issue
%    added ERL limit and Nd
% r227 adding Pmax/Vf and peak of eq pulse, fixed issue with rx testing
%     INC_PACKAGE=0 not fully supported message
%     if N=0 use TDR_duration
%     red display text for fail ERL and COM
% r228 fixed ERL pass fail report, default Grr_limit to 1
% r230 add rx ffe
% r231 change crosstalk noise to icn like to speed things up
% r231 change default OP.impulse_response_truncation_threshold to 1e-3 from
% 1e-5mof-
% r232 fix default for Rx eq so old spead sheets work
% r234 fix inadvertent typo for clause 120e ctle and problem with TXffe loop time reduction
% r235 adding dfe quantization changed to normalized DFE taps reported
% r236 adding ffe gain loop and resample after RxFFE
% r240 added output for C2M and setting defaults for some FFE eq
% r241 force FFE main cursor to 1 and remove sum of taps = 1
% r250 adding more complex package
% r251 post cursor fix for DFE in force() and ffe backoff
% r251 remove TDR threshold noise filter
% r252 add rx FFE filter to receiver noise filter
% r252 change ICN in the xtk noise calculation to end at fb rather than fb/2
% r253 a few bug fixes in force from i indexing and for no ffe postcursors
% r254 precursor check fix in optimize_fom % mod fix in force
% r254 help to align columns in csv file
% r254 accept syntax for 2 tline flex package model
% r256 speed up optimize FOM
% r256 fix problem reading in config file from q/a
% r256 added code from Yasou Hidaka for reading in parameter an and printing out noise
% r257 fixed extrapolation of channel with lower bandwidths in s21_to_impulse_DC
% r257 in get_xtlk_noise in optimize_FOM: reomove crosstalk double counting and apply TXFFE is  FEXT
% r258 EXE_MODE switch 12/21 0:legacy 1:fast 2:superfast
% r258 CDR switch 'MM' or 'mod-MM'
% r258 correction for asymentirc tx/Rx packages
% r258 revamped display results display window
% r259 fix problem if Min_VEO is set in spreadsheet.
% r259 fix problem in optimize_FOM. get_xtlk_noise need to have 3 output
%      parameter else only FEXT is considered for FOM.zhilei huang 01/11/2019
% r259 putting COM_db and IL last in output to terminal
% r259 msgtext change to msg for C2C case other cases not vetted but not presently used
% r259 use N_bx for ERL rather than Nb (ndfe))
% r259 added TDR_W_TXPKG which performs TDR and ERL with the Tx package added
% r260 r259 used rd for the reciever to terminate the package. It was changed to the rd of the transmitter
% r260 used eta_0 PSD equation for sigma_n
% r260 fix IL graph legend to w/pkg and Tr
% r260 define tfx for each port
% r262 fx parameter passing parsing for mod_string revert to 2.57 no COM computational impact
% r262 Report estimate for DER for channel (Yasuo 2/30/19)
% r262 reset on exit default text interpreter to tex
% r262 localize run timer (John Buck 1/17/19)
% r262 set db as internal function in force to avoid tool box
% r262 changed loop for Grr and Gloss in get_tdr so that nbx and tfx works when beta x = 0
% r263 added to output_args RL structure and report "struct" in csv file
% r264 added EW estimate
% r266 using unequalized IR for Vf and Vf to compute ratio of Vp/Vf
% r267 added floating taps with param.N_bf, param.N_bg, param.N_bmax, param.bmaxg.  groups not used for ERL
% r268 added sequential/co-optimization switch for floating tap banks OP.FT_COOP default 0 i.e. sequential
% r269 changed  param.N_bmax to  param.N_f
% r270 implement JingBo Li's and Howard Heck's floating tap method
% r270 modification by Adam Healey for Ls and Cb termination (aka t-coil emulation)
% r270 added c_0 and c_1 for CA in add_brd
% r272 fixed version syntax problem in output_args RL report
% r272 fixed eye width computation problem crosstalk was missed in pervious versions
% r272 removed eye width report if doing a Rx calibration
% r273 better alignment and control for ICN reporting
% r273 fixed PSXTK graph
% r275 fixed delay adjustment for ERL/TDR in get_TDR (Adam Healey 09/06/2019)
% r276 go back to reporting channel IL results (output_args.IL_dB_channel_only_at_Fnq) with board added read_s4p_files (as in r270)
% r276 chdata(i).Aicn=param.a_icn_fext should have been chdata(i).Aicn=param.a_icn_next for the next selection. Since in most spec's they are the same there is little no impact in results
% r276 test for output_args for isfield(chdata(1),'sdd22_raw')
% r276 change divisor for ICN and FOM_ILD to param.f2 from param.fb, may raise ICN and ILD value reported in r275
% r276 C_1 was instantiated as C_0. This was fixed
% r276 fixed rounding problem in reporting of loss at f_nq
% r276 power limit (RSS) for tail DFE taps (B_float_RSS_MAX, N_tail_start)
% r277 added nv for deterining steady state voltage for fitting compatibility
% r278 added b_min to support asymmetric bmax
% r278 added kappa1 and kappa2 to scale package to channel reflection for ERL experiments
% r278 added keyword OP.SHOW_BRD which includes added board in TDR and ERL
% r292 speed up for FOM search (Adee Ran) implemented by Adam Gregory.
% r292 param.LOCAL_SEARCH set to is the heuristic step distance keyword is 'Local Search'
% r292 fixing TDR for different impedance references in get_TDR and s2p file compatibility
% r292 eq. 93A-19 and 93-20 code implementation  bug when include .3by change% to fix edge rate equation 93A-46 (h_T). no effect if Rd=50 or IL > 5 dB
% r292 H_t implemented  in s21_pkg
% r292 plot and report for die to die IL remove the Tr effect "IL with pkgs & Tr filter" goes to "IL with pkgs"
% r292 add GDC_MIN to optimize_FOM
% r293 fix if ndfe-0 and ERL only and s2p issue
% r293a investigate the Tukey filtering
% r293a if fix if bmin is missing
% r294 fix problems reading s2p files for ERL computation
% r294 align Tukey_Window with .3ck definition for ERL and TDR computations
% r294 add parameter param.Noise_Crest_Factor. Default is not to use
% r294 add gdc and gdc2 range limitations
% r295 add VEC Pass threshold
% r295 removed close force all. Tagged all figures with "COM"
% r295 consolidated print in new function "end_display_control"
% r295 report pre/pmax for Txffe
% r295 speed up test cases by not re-reading in s4p files
% r297 add  provisions for AC_CM_RMS for through CM (experimental)
% r299 add keyword T_O (param.T_O) and  samples_for_C2M (parsm.samples_for_C2M) for new C2M VEC and EH computations
% r310 refine VEC and EH for C2M from Adam Gregory in
% r315 added keyword for Bessel_Thomson and Butterworth(default) filter. 
%      cdf_to_ber_contour,COM_eye_width,combine_pdf_same_voltage_axis,
%      optimize_fom_for_C2M. pdf_to_cdf, conv_fct_MeanNotZero, and get_pdf_full
% r311 added RILN
% r314 when T_O is not zero 3 eyes are used to compute VEC and VEO
% r315 Bessel_Thomson keyword is added mostly for measuring  Pmax, Vf, and SNDR
% r316 remove DC computation for RX Calibration loops
% r317 for SAVE_TD to include EQ and unEQ FIR
% r317 clean up bessel thomson and butterworth filter logic for ERL and normal COM 
% r318 if min_VEO_test fails to find a solution the loop is restarted with min_VEO_test to near zero. Makes sure COM returns results 
% r320 fixed RX_CALIBRATION which was broken in r310
% r320 speed up for C2M by moving managing optimize loop distribution of computations
% r320 for C2M added Gaussian window keyword, Gaussian_histogram_window, for T_O and keyword, QL which is  at Q limit at +/-T_O
% r320 removed external feature and replace with TDMODE
% r320 added TDMODE which allows for the use of pulse resonance files (CSV) instead of s4p files
% r330 changed FOM ILN to use a complex fit and compute FOM_ILN in the time domain330 added tfx to N for ERL
% r335 fixed typo in when processing the bessel thompson filter option
% r335 process in CD mode instead of DC mode to get CM noise at Rx
% r335 compute and report CD_CM_RMS
% r335 fixed where output_arg is save i.e. move to end
% r335 refine interp_Sparam to do zero fill instead of extrapolation
% r335 change raw IL plot to not include boards
% r335 set T_0 to zero if not C2M
% r335 change for s parameter interp: check fit sigma, if not OK zero fill 
% r335 added actual sdd12 (instead fo mirroring sdd21) to s21_pkg and 21dc_pkg and read_s4p_files
% r335 TD_RILN changes from Hansel Dsilva
% r335 Fixed sigma_N for RxFFE
% r335 added more to self documenting keyword capability from read_ParamConfigFile and xls_parameter routines
% r335 added c(2) and C(3) back to read_ParamConfigFile
% r335 Optimize_loop_speed_up keyword option added.  Mostly speeds up c2m(vsr)
% r335 corrected GDC_MIN per 0.3ck D2.3
% r335 sigma_r replaces Qr which replaced QL for Gaussian histogram window
% r340 fix for when post cursor taps 2 and 3 are used (from Matt Brown)
% r370 speed up
% r370 fix for floating tap missing locations
% r370 variable Tx FFE taps 
% r370 package die load with ladder circuit
% r370 mods for SNDR_tx exporation using keyword SNR_TXwC0
% r380 fix for Rx Calibaration (error introduced going from 3.4 to 3.7)
% r380 added capabablity to enable a raised cosine Rx filter0
% r380 keyword added: RC_Start, RC_end, Raised_Cosine
% r380 added plot for VTF 
% r385 added capability for additional Tx FFE per package
% r385 keyword added: PKG_Tx_FFE_preset default is 0 i.e. noop
% r385 SAVE_CONFIG2MAT set to 0 as default i.e. don't create a conifig mat file(
% r388 Adjusted Rx caliberation for CL 162 i.e. adding Tx noise (sigma hn) instead of Rx noise line
% r389 Improvement by A. Ran for reporting loss at Nq
% r389 Fixed typo: changed VIM to VMP 
% r400 fixed PR with zero pad extension
% r400 keyword MLSE and SNRADJ_EQUA for future work
% r400 replaced function db with instances of 20*log10(abs(...))
% r410 widen voltage distriution for normal_dist doubled max Q
% r410 improve reading in of config files
% r410 renormalize s-parameter if not 50 ohm ref
% r410 reference for RXFFE changed to MM from UI+zero first precursor
% r410 remove RL from output_args bc not need and too much storage allocation
% r410 s21^2 changed to s12*s21  in s21_pkg. Corrected VTF needed for non-passive sparameters
% r420 updated equalization figures. if Rxffe is use a subplot of Rx FFE taps is graphed
% r420 updated equalization figures. Now separate per pkg case in optimize_fom
% r420 updade force to account for pulse responces with short delays in force
% r420 added Tx/Rx p/n skew with keywords Txpskew, Txnskew, Rxpskew, Pxnskew
% r420 add common mode outputs: VMC_H_mV and SCMR_dB from CDF of CD PR and DD PR
% r420 fixed and added control for RXFFE_TAP_CONSTRAINT and RXFFE_FLOAT_CTL 
% r420 Wiener-Kofp MMSE optimization for RxFFE 
