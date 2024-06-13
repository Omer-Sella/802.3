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