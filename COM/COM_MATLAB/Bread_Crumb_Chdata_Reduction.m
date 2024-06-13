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