% Author: Arjun Chakrawal
% License:This code is released under the CC-By4.0 License.
% Contact:
%   For any questions or concerns, please contact the author at
%   arjun.chakrawal@natgeo.su.se
%%
clear all
close all
clc
fid = fopen('simulation_progress.txt', 'a'); % open text file in append mode
fprintf(fid, 'Starting run_all_fitting.m\n');
fprintf(fid, '\n');%%
search_string = 'fitting';
m_files = dir('*.m');
matching_files = {};
for i = 1:numel(m_files)
    if contains(m_files(i).name, search_string)
        matching_files{end+1,1} = m_files(i).name;
    end
end

matching_files=[matching_files;"Lig_decomposition_starts_time.m"];

for i = 1:length(matching_files)
    time = datestr(clock, 'YYYY/mm/dd HH:MM:SS:FFF');
    fprintf(fid,"Starting "+ matching_files{i}+...
        ".........%23s\n", time);
    run_test(matching_files{i})
    s=findobj('type','legend'); delete(s)
    time = datestr(clock, 'YYYY/mm/dd HH:MM:SS:FFF');
    fprintf(fid,"Finished "+ matching_files{i}+...
        ".........%23s\n", time);
    fprintf(fid, '\n');
end

fprintf(fid, 'All simulations finished!!\n');
fclose(fid); % close the text file

%%
function run_test(fname)
run(fname)
end


