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