
clear all
close all
clc
fid = fopen('simulation_progress.txt', 'a'); % open text file in append mode
fprintf(fid, 'Starting run_all_fitting.m\n');
fprintf(fid, '\n');

search_string = 'fitting';
m_files = dir('*.m');
matching_files = {};
for i = 1:numel(m_files)
   if contains(m_files(i).name, search_string)
      matching_files{end+1,1} = m_files(i).name;
   end
end



for i = 1:length(matching_files)
   time = datestr(clock, 'YYYY/mm/dd HH:MM:SS:FFF');
   fprintf(fid,"Starting "+ matching_files{i}+...
      ".........%23s\n", time);
   run_test(matching_files{i})
   s=findobj('type','legend'); delete(s)
%    run_test('test.m')
   time = datestr(clock, 'YYYY/mm/dd HH:MM:SS:FFF');
   fprintf(fid,"Finished "+ matching_files{i}+...
      ".........%23s\n", time);
   fprintf(fid, '\n');
%    pause
end

fprintf(fid, 'All simulations finished!!\n');
fclose(fid); % close the text file

%%
function run_test(fname)
run(fname)
    % Copy the contents of test.m here, but remove the "clear all" statement
    % and any other unnecessary code that you don't want to run
end