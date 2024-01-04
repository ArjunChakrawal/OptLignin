


% Define the search and replace strings
% search_str1 = 'rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,tau,Ctau];';
% replace_str1 = 'rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];';
search_str2= 'final_table(id, 9:24) = array2table(par_);';
replace_str2= 'final_table(id, 9:26) = array2table(par_);';

% Get a list of all the files in the directory
file_list = dir('*.m');
search_string = 'fitting';
m_files = dir('*.m');
matching_files = {};
for i = 1:numel(m_files)
    if contains(m_files(i).name, search_string)
        matching_files{end+1,1} = m_files(i).name;
    end
end
% Loop over each file and replace the text
for i = 1:length(matching_files)
    file_path = matching_files{i};
    file_contents = fileread(file_path);
%     new_contents = strrep(file_contents, search_str1, replace_str1);
    new_contents = strrep(file_contents, search_str2, replace_str2);
%     new_contents = strrep(new_contents, 'final_table(id, 9:23)', 'final_table(id, 9:24)');
    fid = fopen(file_path, 'w');
    fwrite(fid, new_contents);
    fclose(fid);
end

