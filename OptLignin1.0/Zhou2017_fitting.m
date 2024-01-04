clear all
clc
close all
%%

set_up

disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";


%%
data=readtable('data/Zhou etal 2017.xlsx','Range','A4:G56');
CN0=51.68;
AISC0=0.219;
treatementName = unique(data.treatment);

%%
par_=[];
for i =1:length(treatementName)
    temp=data(data.treatment==treatementName(i),:);
    obs_data=[];
    obs_data.tobs  = temp.day;
    obs_data.Ct_obs = temp.CG;
    obs_data.Co_obs  = aromatic_fraction_inAIS(temp.LigninCG);
    tnorm = obs_data.tobs;
    y=log(obs_data.Ct_obs./obs_data.Ct_obs(1));
    final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
    k=(y(end)-y(1))./tnorm(end);
    terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
    n=terminalTime./obs_data.tobs(end);
    t=0:tnorm(end)*n;
    figure(1);clf
    subplot(121)
    plot(obs_data.tobs,obs_data.Ct_obs,'b-o','DisplayName','Litter C'); hold on
    plot(t,obs_data.Ct_obs(1).*exp(k.*t))
    plot(obs_data.tobs,obs_data.Co_obs,'DisplayName','AromaticC'); hold on
    xlabel("time"); ylabel('gC/bag')
    % title(erase( soilname{i} , "_" ))z
    subplot(122)
    plot(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1),'-o'); hold on
    % title(erase( soilname{i} , "_" ))
    xlabel("litter C (fraction of initial)"); ylabel('lignin (normalized to initial value)')
    param.emax = emax_fun(CN0);
    param.CO_0=obs_data.Co_obs(1);
    param.CT_0=obs_data.Ct_obs(1);
    CO_N= max(obs_data.Co_obs./obs_data.Co_obs(1));
    if(CO_N>1.8)
        mo_f=3;
        disp("mo_f=2")
    else
        mo_f = 1;
    end

    param.lb = [0.0001, param.mo*mo_f, 0.1];
    param.ub = [0.04, param.mo*mo_f, 400];


    fig=figure(2);
    init_guess = [0.006, 0.1,400]; % [vh_max, mo,ro];
    [par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
        init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
    param.rmse=rmse;
    param.rsquare=rsquare;
    makeplot_state_space(sol, '',param,obs_data,g,fig)
    exportgraphics(fig, fig_path_all+"Zhou2017_i"+i+".png",'Resolution',100)
    par_ = [par_;i, par,AISC0,CN0,AISC0*CN0,rmse,rsquare];

end
save_fname = excel_path+ "fitted_par_Zhou2017.txt";
save(save_fname,"par_",'-ascii','-double','-tabs')


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









