clear all
clc
close all
%%

set_up
study ="Tu2011";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
%%

SP=readtable('..\data/Tu et al 2011.xlsx','Range','A10:A33');
time=readtable('..\data/Tu et al 2011.xlsx','Range','B10:B33');
C_remaining= readtable('..\data/Tu et al 2011.xlsx','Range','P10:S33');  %gC litter
lignin_C_remaining= readtable('..\data/Tu et al 2011.xlsx','Range','T10:W33');
init_chem = readtable('..\data/Tu et al 2011.xlsx','Range','A3:J6');  %gC litter
treatementName = C_remaining.Properties.VariableNames;
SPname = unique(SP.Type);
%%
par_=[];
fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];
for i = 1:3
   id= find(strcmpi(SP.Type,SPname{i}));
   for j=1:4
      init_CN=init_chem.initC_N(strcmpi(init_chem.Type,SPname{i}));
      param.emax = emax_fun(init_CN);
      obs_data=[];
      obs_data.tobs  = time.day(id);
      obs_data.Ct_obs = C_remaining(id,j).Variables;
      obs_data.Co_obs  = aromatic_fraction_inAIS(lignin_C_remaining(id,j).Variables);
      
      init_guess = [0.003	0.003 0.2	50];
      [par, sol] = find_parameter_lsq(obs_data, param, ...
         init_guess, @ysim_state_space_lsq)  ;
      [r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);
      
      nexttile
      fit_plotting(sol, obs_data);lh = legend('show');lh.Box = "off";
      title("i=" +i)
      CT=sol.y(1,:);Co=sol.y(2,:);
      
      dcodt = diff(sol.y(2,:))./diff(sol.x);
      neg_idx  = find(dcodt<0, 1);
      if(isempty(neg_idx))
         neg_idx=length(sol.x);
      end
      tau=sol.x(neg_idx);
      CTN=CT./CT(1); Ctau=CTN(neg_idx); 
      
      LC0=init_chem.initL_C(strcmpi(init_chem.Type,SPname{i}));
      par_ = [par_; par,LC0,init_CN,LC0*init_CN,...
         r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];
      j
   end
   i
end
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)





