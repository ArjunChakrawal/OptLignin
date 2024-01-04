clear all
clc
close all

%%
set_up
study = "Berg";

disp("this is the g-function with " +scenario_name+" scenario")
excel_path = "est_params\";
fig_path = "fig\";

% Read_data
initial_fraction_of_C = 0.5;
berg_data = readtable('..\data\Berg and McClaugherty 1989 - Suppl Mat.xlsx');
datacode = unique(berg_data.DatasetCode);
init_guess = [0.003	0.003 0.2	50];

%% Berg fitting
close all
ix = 0;
fig = figure;
fig.Position = [100, 100, 1200, 800];
tiledlayout('flow');
tau = [];
Ctau=[];
% param.lb = [0.0001,0.0001, param.mo, 0.001];
% param.ub = [0.04,0.04, param.mo, 100];
par_ = [];
for i = 1:length(datacode)
   id = find(berg_data.DatasetCode == i);
   dataset = berg_data(id(1):id(end), :);
   spcode = unique(dataset.SpeciesCode);
   for j = 1:length(spcode)
       disp("current simulation is i =" +i+" j=" +j)
      ix = ix + 1;
      idsp = find(dataset.SpeciesCode == spcode(j));
      decdata2 = dataset(idsp(1):idsp(end), :);
      id = decdata2.NMg_g < 0;
      decdata2(id, :) = [];
      total_C = (1-decdata2.MassLoss_*0.01)*initial_fraction_of_C; % gC/g litter
      amount_AUR_C = (decdata2.AISMg_gInitialLitter*0.001*0.5).*(1-decdata2.MassLoss_(1)*0.01);  % gligC/g litter
      aromaticC = aromatic_fraction_inAIS(amount_AUR_C); % true lignin
      
      C_2 = decdata2.gC_gCInitial(2);
      % N_2 = decdata2.gN_gNInitial(2);
      if (C_2 < 0.7)
         tsrt = 2;
         temp = 1;
      else
         tsrt = 1;
         temp = 0;
      end
      obs_data = [];
      Time_days_= decdata2.Time_days_(tsrt:end) - decdata2.Time_days_(tsrt);
      obs_data.tobs = Time_days_;
      obs_data.Ct_obs = total_C(tsrt:end);
      obs_data.Co_obs = aromaticC(tsrt:end);
      
      param.CO_0 = obs_data.Co_obs(1);
      param.CT_0 = obs_data.Ct_obs(1);
      
      AISC0=decdata2.gLigninC_gC(tsrt);
      CN0=1./decdata2.gN_gC(tsrt);
      LN0=decdata2.gligC_gN(tsrt);
      
      param.emax = emax_fun(CN0);
      
      init_guess = [0.003	0.0003 0.1	50];
      [par, sol] = find_parameter_lsq(obs_data, param, ...
         init_guess, @ysim_state_space_lsq)  ;
      [r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);
      
      nexttile
      fit_plotting(sol, obs_data)
      str="i=" +i+", j=" +j;
      str2="(" +decdata2{1, 1}+":" +decdata2{1, 2}+")";
      title({str,str2},'FontSize',8);
      CT=sol.y(1,:);Co=sol.y(2,:);
      
      dcodt = diff(sol.y(2,:))./diff(sol.x);
      neg_idx  = find(dcodt<0, 1);
      if(isempty(neg_idx))
          neg_idx=length(sol.y);
      end
      tau=sol.x(neg_idx);
      CTN=CT./CT(1); Ctau=CTN(neg_idx); 
      par_ = [par_; par, AISC0, CN0, LN0, r2,...
         rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];
%       if(r2_co<0)
%           error("rsquare negative")
%       end
      
   end
end

%%
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)






