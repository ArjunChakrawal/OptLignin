clear all
clc
close all

%
set_up
excel_path = "est_params\" +scenario_name + "\";
fit_results_table = readtable("data\data summary.xlsx");

search_string = '_final';
files = dir(excel_path+'*.xlsx');
matching_files = {};
for i = 1:numel(files)
    if contains(files(i).name, search_string)
        matching_files{end+1,1} = files(i).name;
    end
end

id =strcmp(matching_files,'Preston09NMR_final.xlsx');
matching_files(id)=[];
for i = 1:length(matching_files)
    fname = matching_files{i};
    data = readtable(excel_path+fname);
    sz = size(data, 2);
    %     newstr = extractAfter(fname,"par_");
    studyname = erase(fname, "_final.xlsx");

    id = strcmp(fit_results_table.study_name, studyname);
    if (size(data, 2) == 23)
        fit_results_table(id, 9:29) = data(:, 3:end);
    else
        fit_results_table(id, 9:29) = data(:, 2:end);
    end
    i

    %     for i=1:sz
    %         c=[c;{studyname}];
    %     end
    %     pause
end

aromaticC0=fit_results_table.AURC0;
for i = 1:size(fit_results_table, 1)
    id = fit_results_table.LigninMeasurementMethod_Klasson_1_AcidDetergentFiber_2_CuO_VSC_(i) ~= 3;
    if (id)
        aromaticC0(i) = 0.2 * fit_results_table.AURC0(i);
    end
end

fit_results_table.aromaticC0=aromaticC0;   
fit_results_table.("AC:N0") = fit_results_table.aromaticC0 .* fit_results_table.CN0;

%%
T=fit_results_table;
T.Ntreat = T.TreatmentControl_0_Nitrogen_1_2_3_4;
T.Ntreat(T.Ntreat~=0)=1;

controldata = T(T.Naddexp_1_else_0==0,:);
controldata.logRR_vh=ones(size(controldata,1),1)*-999;
controldata.logRR_avgvo=ones(size(controldata,1),1)*-999;
controldata.logRR_ro=ones(size(controldata,1),1)*-999;
controldata.deltaTau=ones(size(controldata,1),1)*-999;

NAdddata = T(T.Naddexp_1_else_0==1,:);

study = unique(NAdddata.study_name);
logRR_vh= [];logRR_avgvoO= [];logRR_ro= [];deltaTau= [];

for i =1:length(study)
    
    temp = NAdddata(strcmp(NAdddata.study_name,study(i)),:);
    SP = unique(temp.Species);
    for j=1:length(SP)
        id=strcmp(temp.Species,SP(j));
        temp2 = temp(id,:);
        id2=strcmp(NAdddata.Species,SP(j));
        NAdddata.logRR_vh(id2) = log(temp2.vhmax/temp2.vhmax(temp2.Ntreat==0));
        NAdddata.logRR_avgvo(id2) = log(temp2.avg_vo/temp2.avg_vo(temp2.Ntreat==0));
        NAdddata.logRR_ro(id2) = log(temp2.ro/temp2.ro(temp2.Ntreat==0));
        NAdddata.deltaTau(id2) = temp2.LigDecStartDay-temp2.LigDecStartDay(temp2.Ntreat==0);
        disp([i,j]);
    end
end

%%
mydata= [controldata;NAdddata];
delete(excel_path+"fit_results_table.xlsx");
writetable(mydata,excel_path+"fit_results_table.xlsx")


