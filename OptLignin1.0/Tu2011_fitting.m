clear all
clc
close all
%%

set_up
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";

param.lb=[0.0005, param.mo,0.01];
param.ub=[0.04, param.mo,400];
%%

SP=readtable('data/Tu et al 2011.xlsx','Range','A10:A33');
time=readtable('data/Tu et al 2011.xlsx','Range','B10:B33');
C_remaining= readtable('data/Tu et al 2011.xlsx','Range','P10:S33');  %gC litter
lignin_C_remaining= readtable('data/Tu et al 2011.xlsx','Range','T10:W33');
init_chem = readtable('data/Tu et al 2011.xlsx','Range','A3:J6');  %gC litter
treatementName = C_remaining.Properties.VariableNames;
SPname = unique(SP.Type);
par_=[];
for i = 1:3
    id= find(strcmpi(SP.Type,SPname{i}));
    for j=1:4
        init_CN=init_chem.initC_N(strcmpi(init_chem.Type,SPname{i}));
        param.emax = emax_fun(init_CN);
        obs_data=[];
        obs_data.tobs  = time.day(id);
        obs_data.Ct_obs = C_remaining(id,j).Variables;
        obs_data.Co_obs  = aromatic_fraction_inAIS(lignin_C_remaining(id,j).Variables);

        tnorm = obs_data.tobs;
        f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
        final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
        k=f.b ;
        terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
        n=terminalTime./obs_data.tobs(end);
        tt=0:tnorm(end)*n;

        param.CO_0=obs_data.Co_obs(1);
        param.CT_0=obs_data.Ct_obs(1);

        figure(1);
        subplot(121)
        plot(obs_data.tobs,obs_data.Ct_obs,'o-','DisplayName','TC'); hold on
        plot(tt,obs_data.Ct_obs(1).*exp(k.*tt))
        scatter(terminalTime,obs_data.Ct_obs(1).*exp(k.*terminalTime),100,'Marker','*', ...
            'MarkerEdgeColor','red')
        plot(obs_data.tobs,obs_data.Co_obs,'o-','DisplayName','ligninC'); hold on
        xlabel("time"); ylabel('gC/bag')
        subplot(122)
        scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
        xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
        suptitle(SPname{i}+" "+ treatementName{j})


        
%         init_guess = [0.0036, 0.0001,0.01]; % [vh_max, mo,ro];
%         [ocp,sol] =  opt_con(param,g,init_guess,obs_data.tobs(end)*n);
%         makeplot_state_space(sol, '',param,obs_data,g,figure)

        fig=figure(2);
        init_guess = [-k, 0.05,2]; % [vh_max, mo,ro];
        [par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
            init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
        param.rmse=rmse;param.rsquare=rsquare;
        makeplot_state_space(sol, '',param,obs_data,g,fig)
        exportgraphics(fig, fig_path_all+"Tu2011_"+SPname{i}+"_"+ ...
            treatementName{j}+".png",'Resolution',100)

        LC0=init_chem.initL_C(strcmpi(init_chem.Type,SPname{i}));
        par_ = [par_;i, j, par,LC0,init_CN,LC0*init_CN,rmse,rsquare];

    end

end
save_fname = excel_path+ "fitted_par_Tu2011.txt";
save(save_fname,"par_",'-ascii','-double','-tabs')

