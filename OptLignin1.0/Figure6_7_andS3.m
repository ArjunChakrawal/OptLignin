
clear all;
close all;clc
set_up
excel_path = "est_params\"+scenario_name+"\";
fname= excel_path+"fit_results_table.xlsx";
rawdata=readtable(fname);

rawdata = renamevars(rawdata, {'LigninMeasurementMethod_Klasson_1_AcidDetergentFiber_2_CuO_VSC_',...
    'LitterTypeGrass_1Leaf_2Needle_3Roots_4Wood_5Lichen_6','LigDecStartDay',...
    'ClimateTundra_1Boreal_2CoolTemperate_3WarmTemperate_4Tropical_5',...
    'vo_at_Ligdecomposition_Starts','AC_N0','aromaticC0'},...
    {'Lig_method','LitterType','tau','climate','voAtTau','AN0','ARC0'});

id=rawdata.r2_co<-0.25;
sum(id)
rawdata(id,:)=[];
study =  unique(rawdata.study_name);
id = 1:length(study);

for i =1:length(study)
    idd=strcmp(rawdata.study_name,study{i});
    rawdata.studyID(idd)=id(i);
end

c= min(rawdata.tau(rawdata.tau>0))/2;
rawdata.tautf = log(rawdata.tau+c);

rawdata.MAT_S = normalize(rawdata.MATC);
rawdata.MAP_S = normalize(rawdata.MAPMm);

rawdata.logro=log(rawdata.ro);
rawdata.logvh=log(rawdata.vhmax);
rawdata.logavgvo=log(rawdata.avg_vo);
rawdata.logmax_vo=log(rawdata.max_vo);
rawdata.logvoAtTau=log(rawdata.voAtTau);

rawdata.logARC0=log(rawdata.ARC0);
rawdata.logCN0=log(rawdata.CN0);
rawdata.logAN0=log(rawdata.AN0);

rawdata.CN0_S=normalize(rawdata.CN0);
rawdata.AN0_S=normalize(rawdata.AN0);
rawdata.ARC0_S = normalize(rawdata.ARC0);

R = 8.314e-3; %kJ mol-1 K-1
col={'studyID','LitterType','Ntreat','MATC','MAPMm','CN0','ARC0','AN0','tau','vhmax','avg_vo','ro','max_vo','voAtTau'};
tempdata=rawdata(:,col);
delete("data_forRStudio.xlsx");
writetable(tempdata,"data_forRStudio.xlsx")

col={'studyID','MAT_S','MAP_S','CN0_S','ARC0_S','tautf','logro','logvh','logavgvo','logmax_vo'};

%% final regression
Rpath = 'C:\PROGRA~1\R\R-43~1.2\bin\x64';
RscriptFileName = 'rsquared_values.R';
RunRcode(RscriptFileName, Rpath)

rawdata.studyID = categorical(rawdata.studyID);
rawdata.LitterType = categorical(rawdata.LitterType);

target ={'tautf','logro','logvh','logavgvo','logmax_vo'};
pred_formula='~ MAT_S*ARC0_S+MAT_S*CN0_S+CN0_S*ARC0_S  + (1|studyID)';
[lme,~] =checkAssumptions(rawdata,target,pred_formula,'lme','FitMethod','ML');
yticklabels({'log($\tau$+K)','log({$r_O$})','log({$v_H$})', ...
    'log($\bar{v}_O$)','log(max(${v}_O$))'})
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
xticklabels({'MAT_S','\it{CN_{0,S}}','{\itARC_{0,S}}','MAT_S\times{\itCN_{0,S}}' ...
    ,'MAT_S\times{\itARC_0_S}','\itCN_{0,S}\times{\itARC_{0,S}}',...
    "({\itr^2_{marg} , r^2_{cond}})"})

title('')
exportgraphics(gcf, 'results\Figure6.png', Resolution=300)

% Temperature sensitivity Q10

beta_vh =lme{3}.Coefficients.Estimate(2)/std(rawdata.MATC);
beta_avgvo =lme{4}.Coefficients.Estimate(2)/std(rawdata.MATC);
beta_maxvo =lme{5}.Coefficients.Estimate(2)/std(rawdata.MATC);


Q10_vh = exp(10*beta_vh);
Q10_avgvo = exp(10*beta_avgvo);
Q10_maxvo= exp(10*beta_maxvo);
coeff_T_Allison2018_hydro = log(mean([1.6,2.25]))*std(rawdata.MATC)/10;
coeff_T_Allison2018_ligno = log(1.4)*std(rawdata.MATC)/10;

%% interaction plots
close all
cmap = [255, 102, 102; % red
    0, 153, 153]/255; % blue
numLevels = 5;
LC = interp1(linspace(-1, 1, 2), cmap, linspace(-1, 1, numLevels));
fig=figure(Color='w');fig.Position=[200 200 710 458 ];
tiledlayout('flow', "TileSpacing","compact"); 

sMAT = linspace(min(rawdata.MAT_S),max(rawdata.MAT_S),50)';
sCN0=0;
sARC0 = [-2,-1,0,1,2]';

mod = lme{1};
for i = 1:length(sCN0)
    nexttile
    for j =1:length(sARC0)
        X = ones(length(sMAT),7);
        X(:,2)=sMAT;X(:,3)=sCN0(i);X(:,4)=sARC0(j);
        X(:,5)=sCN0(i).*sMAT;X(:,6)= sARC0(j).*sMAT;X(:,7)=sARC0(j).*sCN0(i);
        Z= X*mod.Coefficients.Estimate;
        plot(sMAT, Z, Color=LC(j,:), LineWidth=2,DisplayName=""+sARC0(j)); hold on 
    end
    title("{\it CN_{0,S}} = "+sCN0(i))
    xlabel("MAT_S");
    ylabel("log ($\tau$+K)", Interpreter='latex');
    axis tight;
    ylim([-1,10]);
    text(0.05,0.95,"(a)", 'Units','normalized')
end

mod = lme{2};
for i = 1:length(sCN0)
    nexttile
    for j =1:length(sARC0)
        X = ones(length(sMAT),7);
        X(:,2)=sMAT;X(:,3)=sCN0(i);X(:,4)=sARC0(j);
        X(:,5)=sCN0(i).*sMAT;X(:,6)= sARC0(j).*sMAT;X(:,7)=sARC0(j).*sCN0(i);
        Z= X*mod.Coefficients.Estimate;
        plot(sMAT, Z, Color=LC(j,:), LineWidth=2,DisplayName=""+sARC0(j)); hold on 
    end
    title("{\it CN_{0,S}} = "+sCN0(i))
    xlabel("MAT_S");
    ylabel("log ($r_O$)", Interpreter='latex');
    axis tight;ylim([-1,4]);
    text(0.05,0.95,"(b)", 'Units','normalized')
end
lh=legend('show'); lh.Box="off";lh.Location='bestoutside';
title(lh,'{\it ARC_{0,S}}')

mod = lme{3};
sCN0 =[-0.5,0,0.5,1,2];
sARC0 = 0;


for i = 1:length(sARC0)
    nexttile
    for j =1:length(sCN0)
        X = ones(length(sMAT),7);
        X(:,2)=sMAT;X(:,3)=sCN0(j);X(:,4)=sARC0(i);
        X(:,5)=sCN0(j).*sMAT;X(:,6)= sARC0(i).*sMAT;X(:,7)=sARC0(i).*sCN0(j);
        Z= X*mod.Coefficients.Estimate;
        plot(sMAT, Z, Color=LC(j,:), LineWidth=2,DisplayName=""+sCN0(j)); hold on 
        axis tight;ylim([-2.5,1]);
    end
    title("{\it ARC_{0,S}} = "+sARC0(i))
    xlabel("MAT_S");
    ylabel("log $({v}_H)$", Interpreter='latex');
    text(0.05,0.95,"(c)", 'Units','normalized')
end
plot(sMAT,sMAT*coeff_T_Allison2018_hydro,'--',LineWidth=2,Color='k',...
    DisplayName='Allison et al. 2018')

mod = lme{4};p=[];
for i = 1:length(sARC0)
    nexttile
    for j =1:length(sCN0)
        X = ones(length(sMAT),7);
        X(:,2)=sMAT;X(:,3)=sCN0(j);X(:,4)=sARC0(i);
        X(:,5)=sCN0(j).*sMAT;X(:,6)= sARC0(i).*sMAT;X(:,7)=sARC0(i).*sCN0(j);
        Z= X*mod.Coefficients.Estimate;
        p(j)=plot(sMAT, Z, Color=LC(j,:), LineWidth=2,DisplayName=""+sCN0(j)); hold on 
        
    end
    title("{\it ARC_{0,S}} = "+sARC0(i))
    xlabel("MAT_S");
    ylabel("log $(\bar{v}_O)$", Interpreter='latex');
%     ylabel("log (max(${v}_O$))", Interpreter='latex');
    axis tight;ylim([-2.5,1]);
    text(0.05,0.95,"(d)", 'Units','normalized')
end
p(6)=plot(sMAT,sMAT*coeff_T_Allison2018_ligno-0.5,'--',LineWidth=2,Color='k',...
    DisplayName='Allison et al. 2018');

lh=legend('show',p); lh.Box="off";lh.Location='bestoutside';
title(lh,'{\it CN_{0,S}}')
exportgraphics(fig, 'results\Figure7.png', Resolution=300)
%% Temperature sensitivity Activation Energy 
close all
rawdata.Tinv = 1./(rawdata.MATC+273);
rawdata.TinvS = normalize(rawdata.Tinv );
pred_formula='~ TinvS*ARC0_S+ TinvS*CN0_S+ CN0_S*ARC0_S  + (1|studyID)';
target ={'logvh','logavgvo','logmax_vo'};
[lme,tbl] =checkAssumptions(rawdata,target,pred_formula,'lme');

Ea_vh=-lme{1}.Coefficients.Estimate(4)/std(rawdata.Tinv)*R; %kJ mol-1
Ea_avgvo=-lme{2}.Coefficients.Estimate(4)/std(rawdata.Tinv)*R; %kJ mol-1
Ea_maxvo=-lme{3}.Coefficients.Estimate(4)/std(rawdata.Tinv)*R; %kJ mol-1

Ea_vh_se=lme{1}.Coefficients.SE(4)/std(rawdata.Tinv)*R; %kJ mol-1
Ea_avgvo_se=lme{2}.Coefficients.SE(4)/std(rawdata.Tinv)*R; %kJ mol-1
Ea_maxvo_se=lme{3}.Coefficients.SE(4)/std(rawdata.Tinv)*R; %kJ mol-1


%% Figure S3 pairwise scatter plot 
close all
LC=summer(3);
rawdata.MAPMm=rawdata.MAPMm*0.001;
% rawdata.tau=rawdata.tau./365;
newdata = renamevars(rawdata, {'MATC','MAPMm','logCN0','ARC0','logAN0',...
    'tau','logvh','logavgvo','logmax_vo','logro'},...
    {'MAT [$^o$C]','MAP[m]','log(CN$_{0})$','ARC$_{0}$','log(AN$_{0})$',...
    '$\tau$','log ($v_H$)','log $(\bar{v}_O)$', ...
    'log(max(${v}_O))$','log ({$r_O$})'});

cols_to_plot= {'MAT [$^o$C]','MAP[m]','log(CN$_{0})$','ARC$_{0}$','log(AN$_{0})$',...
    '$\tau$','log ($v_H$)','log $(\bar{v}_O)$', ...
    'log(max(${v}_O))$','log ({$r_O$})'};
scatterplot_with_correlation(newdata,cols_to_plot,cols_to_plot, 0.05,'',...
    20,'MarkerFaceColor',LC(2,:),'MarkerEdgeColor','none')
fig=gcf;
fig.Position(3:4) = [925,740];
exportgraphics(gcf,"results/FigureS3.png",Resolution=400)

