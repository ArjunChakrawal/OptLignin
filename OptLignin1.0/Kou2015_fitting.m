clear all
clc
close all
%%

set_up

disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";

param.lb=[0.0005, param.mo,1];
param.ub=[0.04, param.mo,400];
%%
data=readtable('data/Kou2015.xlsx','Range','A4:M58');
Littertype = unique(data.type);
treatementName = unique(data.treatment);
par_=[];
for i=1:length(Littertype)
    temp=data(data.type==Littertype(i),:);
    for j =1:length(treatementName)
        temp2=temp(temp.treatment==treatementName(j),:);
        obs_data=[];
        obs_data.tobs  = temp2.day;
        obs_data.Ct_obs = temp2.CInG;
        obs_data.Co_obs  = aromatic_fraction_inAIS(temp2.LigninCInG);
        CN0 = temp2.CN0(1);
        tnorm = obs_data.tobs;
        final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
        f = fit(obs_data.tobs ,obs_data.Ct_obs./obs_data.Ct_obs(1),'exp1');
        k=f.b ;

        terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
        n=terminalTime./obs_data.tobs(end);
        t=0:tnorm(end)*n;
        figure(1);clf
        subplot(121)
        plot(obs_data.tobs,obs_data.Ct_obs,'b-o','DisplayName','Litter C'); hold on
        plot(t,obs_data.Ct_obs(1).*exp(k.*t))
        plot(obs_data.tobs,obs_data.Co_obs,'DisplayName','AromaticC'); hold on
        xlabel("time"); ylabel('gC/bag')
        % title(erase( soilname{i} , "_" ))
        subplot(122)
        plot(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1),'-o'); hold on
        % title(erase( soilname{i} , "_" ))
        xlabel("litter C (fraction of initial)"); ylabel('lignin (normalized to initial value)')

        param.emax = emax_fun(CN0);
        param.CO_0=obs_data.Co_obs(1);
        param.CT_0=obs_data.Ct_obs(1);

        fig=figure(2);
        init_guess = [-k*0.5, param.mo,100]; % [vh_max, mo,ro];
        [par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
            init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
        [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
        param.rsquare = r2;
        param.rmse = rmse;
        makeplot_state_space(sol, '',param,obs_data,g,fig)
        exportgraphics(fig, fig_path_all+"Kou2015_i"+i+"_j"+j+".png",'Resolution',100)
        AISC0 = temp2.L_C(1);
        par_ = [par_;i,j, par,AISC0,CN0,AISC0*CN0,rmse,rsquare];

    end
end
save_fname = excel_path+ "fitted_par_Kou2015.txt";
save(save_fname,"par_",'-ascii','-double','-tabs')


