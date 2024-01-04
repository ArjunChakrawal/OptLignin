clear all
clc
close all
%%

set_up
study ="Zhou2017";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
param.lb = [0.0005,0.0001, param.mo, 0.01];
param.ub = [0.04,0.04, param.mo, 400];

%%
data=readtable('..\data/Zhou etal 2017.xlsx','Range','A4:G56');
CN0=51.68;
AISC0=0.219;
treatementName = unique(data.treatment);

%%
par_=[];
fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];
for i =1:length(treatementName)
    temp=data(data.treatment==treatementName(i),:);
    obs_data=[];
    obs_data.tobs  = temp.day;
    obs_data.Ct_obs = temp.CG;
    obs_data.Co_obs  = aromatic_fraction_inAIS(temp.LigninCG);
    init_guess = [0.003	0.003 0.2	50];
    param.emax = emax_fun(CN0);
    CO_N= max(obs_data.Co_obs./obs_data.Co_obs(1));

    if(CO_N>1.8)
        mo_f=3;
        disp("mo_f=2")
    else
        mo_f = 1;
    end
        param.lb = [0.0001,0.0001, param.mo*mo_f, 0.1];
        param.ub = [0.04,0.04, param.mo*mo_f, 400];

    [par, sol] = find_parameter_lsq(obs_data, param, ...
        init_guess, @ysim_state_space_lsq)  ;
    [r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);

    nexttile
    fit_plotting(sol, obs_data);lh = legend('show');lh.Box = "off";
    title("i=" +i)
    CT=sol.y(1,:);Co=sol.y(2,:);

    dcodt = diff(sol.y(2,:))./diff(sol.x);
    neg_idx  = find(dcodt<0, 1);
    tau=sol.x(neg_idx);
    CTN=CT./CT(1); Ctau=CTN(neg_idx);
    par_ = [par_; par,AISC0,CN0,AISC0*CN0,...
        r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];

end

final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)


%
% n=2
% init_guess = [-k*1.5, 0.4,40]; % [vh_max, mo,ro];
% y1norm=obs_data.Ct_obs./obs_data.Ct_obs(1);
% y2norm= obs_data.Co_obs./obs_data.Co_obs(1);
% ydata = [y1norm;y2norm];
% [ocp,sol] =  opt_con(param,g,init_guess,obs_data.tobs(end)*n);
% ysim=ysim_state_space(init_guess, ocp, obs_data, [param.CT_0, param.CO_0]);
% rmse = sqrt(mean(ysim-ydata).^2);
% SSR =sum((ysim-ydata).^2);
% SST = sum((mean([y1norm; y2norm])-[y1norm; y2norm]).^2);
% rsquare= 1-SSR/SST;
% param.rsquare=rsquare;
% param.rmse=rmse;
% makeplot_state_space(sol, '',param,obs_data,g,fig)



% exportgraphics(fig, fig_path_all+"Hall2020"+soilname{i}+".png",'Resolution',100)









