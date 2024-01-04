clear all;
close all;clc
set_up
excel_path = "est_params\"+scenario_name+"\";
fname= excel_path+"fit_results_table.xlsx";
data=readtable(fname);
%%
Study =  unique(data.Study);
CN0=[];MAT=[];MAP=[];AURC=[];vh=[];avgvo=[];ro=[];tau=[];lignin_method=[];
cliamte=[];
for i =1:length(Study)
    idd=strcmp(data.Study,Study{i});
    data2=data(idd,:);
    if(min(data2.MATC)==max(data2.MATC))
        MAT= [MAT;num2str(min(data2.MATC),'%1.0f')];
        MAP= [MAP;num2str(min(data2.MAPMm),'%1.0f')];
    else
        MAT= [MAT;strcat(num2str(min(data2.MATC),'%1.0f'),"-",num2str(max(data2.MATC),'%1.0f'))];
        MAP= [MAP;strcat(num2str(min(data2.MAPMm),'%1.0f'),"-",num2str(max(data2.MAPMm),'%1.0f'))];
    end
    CN0= [CN0;strcat(num2str(min(data2.CN0),'%1.0f'),"-",num2str(max(data2.CN0),'%1.0f'))];
    AURC= [AURC;strcat(num2str(min(data2.AURC0),'%1.2f'),"-",num2str(max(data2.AURC0),'%1.2f'))];
    vh= [vh;strcat(num2str(min(data2.vhmax),'%1.1E'),"-",num2str(max(data2.vhmax),'%1.1E'))];
    avgvo= [avgvo;strcat(num2str(min(data2.avg_vo),'%1.1E'),"-",num2str(max(data2.avg_vo),'%1.1E'))];
    ro= [ro;strcat(num2str(min(data2.ro),'%1.1f'),"-",num2str(max(data2.ro),'%1.1f'))];
    tau= [tau;strcat(num2str(min(data2.LigDecStartDay),'%1.0f'),"-",num2str(max(data2.LigDecStartDay),'%1.0f'))];

    if(data2.LigninMeasurementMethod_Klasson_1_AcidDetergentFiber_2_CuO_VSC_(1)==1)
        lignin_method=[lignin_method;"Klason"];
    elseif(data2.LigninMeasurementMethod_Klasson_1_AcidDetergentFiber_2_CuO_VSC_(1)==2)
        lignin_method=[lignin_method;"Acid detergent fiber"];
    elseif(data2.LigninMeasurementMethod_Klasson_1_AcidDetergentFiber_2_CuO_VSC_(1)==3)
        lignin_method=[lignin_method;"CuO oxidation"];
    else
        lignin_method=[lignin_method;"Near infrared spectroscopy"];
    end


end

summary_data= table(Study, MAT,MAP,CN0,AURC,vh,avgvo,ro,tau,lignin_method);

writetable(summary_data,"results\data_summary_statistics.xlsx")
