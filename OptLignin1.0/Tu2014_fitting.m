clear all
clc
close all
%%

set_up

disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";

param.lb=[0.0005, param.mo,0.5];
param.ub=[0.04, param.mo,400];
%%

data=readtable('data/Tu_et_al2014.xlsx','Range','B16:O360');
CN0=table2array(readtable('data/Tu_et_al2014.xlsx','Range','K2:K11'));
AISC0=table2array(readtable('data/Tu_et_al2014.xlsx','Range','N2:N11'));
SP=readtable('data/Tu_et_al2014.xlsx','Range','A1:A11');
figure;tiledlayout("flow")
for i =1:10
    nexttile
    tempdata1=data(data.SpeciesCode == i,:);

    %     gscatter(tempdata1.C_C_0_, tempdata1.lignin_lignin_0_, tempdata1.treatment)
    %         hold on; plot(0:1,0:1)
    gscatter(tempdata1.days,tempdata1.C_C_0_,tempdata1.treatment); hold on
    ylim([0,1])
    title(SP.LitterSubstrate{i})
end
%%

treatementName = unique(data.treatment);
SPname = unique(data.SpeciesCode);
par_=[];
for i =1:length(SPname)
    % for i =[6,8]
    tempdata1=data(data.SpeciesCode == SPname(i),:);
    for j =1:length(treatementName)

        tempdata2=tempdata1(tempdata1.treatment==treatementName(j),:);
        C_2 = tempdata2.C_C_0_(2);
        if ( C_2 < 0.7)
            tsrt = 2;
            temp = 1;
        else
            tsrt = 1;
            temp = 0;
        end
        obs_data=[];
        obs_data.tobs  = tempdata2.days(tsrt:end)-tempdata2.days(tsrt);
        obs_data.Ct_obs = tempdata2.CInG(tsrt:end);
        obs_data.Co_obs  = aromatic_fraction_inAIS(tempdata2.ligninCInG(tsrt:end));

        figure(1);clf
        scatter(obs_data.tobs,obs_data.Ct_obs./obs_data.Ct_obs(1));hold on
        scatter(obs_data.tobs, obs_data.Co_obs./obs_data.Co_obs(1));

        CT_N = obs_data.Ct_obs./obs_data.Ct_obs(1);
        iddiff = [false;diff(CT_N)./CT_N(1:end-1) >0.1];
        obs_data.tobs(iddiff)=[];obs_data.Ct_obs(iddiff)=[];obs_data.Co_obs(iddiff)=[];

        plot(obs_data.tobs,obs_data.Ct_obs./obs_data.Ct_obs(1),'o-');hold on
        plot(obs_data.tobs, obs_data.Co_obs./obs_data.Co_obs(1),'o-');
        ylim([0,inf])

        tnorm = obs_data.tobs;
        f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
        final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
        k=f.b ;
        terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
        n=terminalTime./obs_data.tobs(end);
        tt=0:tnorm(end)*n;

        figure(1);
        subplot(121)
        plot(obs_data.tobs,obs_data.Ct_obs,'o-','DisplayName','TC'); hold on
        plot(tt,obs_data.Ct_obs(1)*exp(k.*tt))
        plot(obs_data.tobs,obs_data.Co_obs,'o-','DisplayName','ligninC'); hold on
        xlabel("time"); ylabel('gC/bag')
        subplot(122)
        scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
        xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')

        param.CO_0=obs_data.Co_obs(1);
        param.CT_0=obs_data.Ct_obs(1);
        param.emax = emax_fun(CN0(i));
        init_guess = [-k, 0.1,0.1]; % [vh_max, mo,ro];

        fig=figure(2);clf
        [par,sol,~,~] = find_parameter(obs_data,...
            param,init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
        [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);

        param.rsquare = r2;
        param.rmse = rmse;
        makeplot_state_space(sol, '',param,obs_data,g,fig)

        exportgraphics(fig, fig_path_all+"Tu2014_i"+i+"_j"+j+".png",'Resolution',100)
        par_ = [par_;i,j, par,AISC0(i),CN0(i),AISC0(i)*CN0(i),rmse,r2];
        disp("current simulation is i ="+ i + " j="+ j )
% pause
    end

end
% fittedparTu1(21:24,:)=par_(1:4,:);
% fittedparTu1(29:32,:)=par_(5:8,:);
% fittedparTu2014(34,:)=par_;
%%
save_fname = excel_path+ "fitted_par_Tu2014.txt";
save(save_fname,"par_",'-ascii','-double','-tabs')

% save(save_fname,"fittedparTu2014",'-ascii','-double','-tabs')

