clear all
clc
close all
%%

set_up
study ="Kou2015";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
% param.lb = [0.0005,0.0001, param.mo, 0.01];
% param.ub = [0.04,0.04, param.mo, 400];

%%
data=readtable('..\data/Kou2015.xlsx','Range','A4:M58');
Littertype = unique(data.type);
treatementName = unique(data.treatment);
par_=[];
fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];
for i=1:length(Littertype)
   temp=data(data.type==Littertype(i),:);
   for j =1:length(treatementName)
      temp2=temp(temp.treatment==treatementName(j),:);
      obs_data=[];
      obs_data.tobs  = temp2.day;
      obs_data.Ct_obs = temp2.CInG;
      obs_data.Co_obs  = aromatic_fraction_inAIS(temp2.LigninCInG);
      
      CN0 = temp2.CN0(1);
      param.emax = emax_fun(CN0);
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
      AISC0 = temp2.L_C(1);
      par_ = [par_; [par,AISC0,CN0,AISC0*CN0,...
         r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau]];

   end
end
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)



