clear all;
close all;clc
set_up
excel_path = "est_params\";
fname= excel_path+"data summary_LSQ.xlsx";
%%
data=readtable(fname);

for i = 1:size(data, 1)
    id = data.LigninMeasurementMethod_Klasson_1_AcidDetergentFiber_2_CuO_VSC_(i) == 1 || ...
        data.LigninMeasurementMethod_Klasson_1_AcidDetergentFiber_2_CuO_VSC_(i) == 2 ||...
        data.LigninMeasurementMethod_Klasson_1_AcidDetergentFiber_2_CuO_VSC_(i) == 4;
    if (id)
        data.AISC0(i) = 0.2 * data.AISC0(i);
    end
end
data = renamevars(data, {'AISC0'}, {'aromaticC0'});
data.LN0 = data.aromaticC0 .* data.CN0;
data = renamevars(data, {'LN0'}, {'AC:N0'});


data = renamevars(data, {'LigninMeasurementMethod_Klasson_1_AcidDetergentFiber_2_CuO_VSC_',...
    'LitterTypeGrass_1Leaf_2Needle_3Roots_4Wood_5Lichen_6','LigDecStartDay',...
    'ClimateTundra_1Boreal_2CoolTemperate_3WarmTemperate_4Tropical_5','AC:N0','aromaticC0'},...
    {'Lig_method','LitterType','tau','climate','AN0','ARC0'});
id=data.r2_co<0;
data(id,:)=[];

%%
data.vhmax=data.vhmax*365; % per year
data.vo = data.vo*365; % per year
data.tau = (data.tau)/365; % in year

data.MAT_S = normalize(data.MATC);
data.logro=log10(data.ro);
data.logvh=log10(data.vhmax);
data.logvo=log10(data.vo);

data.logARC0=log10(data.ARC0);
data.logCN0=log10(data.CN0);
data.logAN0=log10(data.AN0);

% data.Ntreat = categorical(data.Ntreat);
data.climate = categorical(data.climate);
data.Study = categorical(data.Study);
data.LitterType = categorical(data.LitterType);
LC=summer(3);

data.Tinv = 1./(data.MATC+273);
data.LNvh=log(data.vhmax);
data.LNavgvo=log(data.vo);
R = 8.314e-3; %kJ mol-1 K-1
%% Welch's t-test https://en.wikipedia.org/wiki/Welch%27s_t-test

lm1=checkAssumptions(data,"LNvh","~Tinv+(1|LitterType)",@fitlme,'StartMethod','random');
slope_lm1= lm1{1}.Coefficients.Estimate(2);
sx_lm1=lm1{1}.Coefficients.SE(2);

ci = coefCI(lm1{1});
EabyR=ci(2,:);
Eavh = EabyR*R; %kJ mol-1
% figure;
% scatter(data.MATC+273, data.LNvh); hold on
% figure;
% scatter(data.Tinv, data.LNvh); hold on
% scatter(data.Tinv, fitted(lme{1}))

lm2=checkAssumptions(data,"LNavgvo","~Tinv+(1|LitterType)",@fitlme,'StartMethod','random');
slope_lm2= lm2{1}.Coefficients.Estimate(2);
sx_lme=lm2{1}.Coefficients.SE(2);

ci = coefCI(lm2{1});
EabyR=ci(2,:);
Eavo = EabyR*R; %kJ mol-1


sdelX = sqrt(sx_lm1^2+sx_lme^2);
test_statistic= (slope_lm1-slope_lm2)/sdelX;
df = sdelX^4/(sx_lm1^4/204 + sx_lme^4/204);
 % two sided for checking hte if the two estimates are different
p_value = 2 * (1 - tcdf(abs(test_statistic), df));

critical_value = tinv(1 - 0.05, df);
 % one sided for checking hte if one slope is higher than the other
p_value = 1 - tcdf(test_statistic, df);
% Compare with p-value
if (p_value < 0.05)
    disp('Reject the null hypothesis. One slope is higher than the other.');
else
    disp('Fail to reject the null hypothesis. The slopes are not significantly different.');
end

%% 
% col={'MAT_S','CN0_S','AN0_S','ARC0_S','tau','ro','logvh','logvo','Ntreat'};
% savname ='rescaled_scatter_poster.png';
% scatterplot_with_correlation(data,col,col, 0.05,savname,...
%     20,'MarkerFaceColor',LC(2,:),'MarkerEdgeColor','none')

newdata = renamevars(data, {'MATC','ARC0','logCN0','logAN0', ...
    'tau','logvh','logvo','logro','logARC0'},...
    {'MAT [^oC]','{\itARC_0} [-]','log ({\itCN_0})','log ({\itARC:N_0})', ...
    '\tau [Y]','log ({\itv_H})','log ({\itavg. v_O})', ...
    'log ({\itr_O})','log ({\itARC_0})'});

LC=summer(3);
savname='All_scatter_plot.png';
cols_to_plot= {'MAT_S','log ({\itARC_0})','log ({\itCN_0})','log ({\itARC:N_0})',...
    '\tau [Y]','log ({\itv_H})','log ({\itavg. v_O})','log ({\itr_O})'};

scatterplot_with_correlation(newdata,cols_to_plot,cols_to_plot, 0.05,savname,...
    20,'MarkerFaceColor',LC(2,:),'MarkerEdgeColor','none')

%% Boxplot

close all;
col={'MAT_S','logCN0','logARC0','logAN0','tau','logvh','logvo','logro'};
T= data(:,col);
T2 = stack(T, col, 'NewDataVariableName', 'Value');
figure;boxplot(T2.Value,T2.Value_Indicator)
%% 
% close all;
clc
target ={'tau','logvh','logvo','logro'};
% pred_formula="~ 1 + MAT_S*logARC0+ MAT_S*logCN0 +(1|study_name)";
pred_formula="~ 1 + MAT_S*logCN0 + logARC0 +(1|study_name) ";

lm2=checkAssumptions(data,target,pred_formula,@fitlme);
title("Linear Regression WITH random effects", ...
    'fontweight','normal')

% close all;clc
% target ={'tau','logro'};
% pred_formula="~ 1 + MAT_S+logCN0 +MAT_S*logvo +(1|study_name)";
% lme=checkAssumptions(data,target,pred_formula,@fitlme,'StartMethod','random');
% title("Linear Regression WITH random effects", ...
%     'fontweight','normal')
% yticklabels({'\tau [year]','log ({\itv_H})','log ({\itavg. v_O})', ...
%     'log ({\itr_O})'})
% xticklabels({'MAT_{s}','log {\itARC_0}','log ({\itCN_0})',...
%     'MAT_{s}*log {\itARC_0}','MAT_{s}*log {\itCN_0}','log {\itCN_0}*log {\itARC_0}'})
exportgraphics(gcf,"results/"+"p_value_final.png",Resolution=300)
