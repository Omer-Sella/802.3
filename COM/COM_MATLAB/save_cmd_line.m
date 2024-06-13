function [ cmd_str ] = save_cmd_line( config_file, chdata, num_fext,num_next, cli_name )
% save commmend string
% for saving from interactive queries


cmd_str=[ cli_name '('  '''' config_file '''' ',' num2str(num_fext) ', ' num2str(num_next) ',' '''' chdata(1).filename ''''];
for i=1:num_next+num_fext
    cmd_str= [cmd_str  ',' '''' chdata(i+1).filename ''''];
end
cmd_str= [ cmd_str ')'];
