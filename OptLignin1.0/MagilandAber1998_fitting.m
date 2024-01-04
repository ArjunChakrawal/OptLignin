clear all
clc
close all
%%

set_up
study="MagilandAber1998";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";

param.lb=[0.0005, param.mo,0.5];
param.ub=[0.04, param.mo,400];
%%
data=readtable('data/Magil and Aber 1998.xlsx','Range','A8:Q242');
%%
figure;tiledlayout("flow")

for i= 1:2
    sitedata = data(data.SiteCode == i,:);
    spcode = unique(sitedata.SpeciesCode);
    for j =1:length(spcode)
        spdata=sitedata(sitedata.SpeciesCode == spcode(j),:);
        treatmentcode = unique(spdata.treatmentCtrl_0_LN_1_HN_2);
        for k =1:length(treatmentcode)
            litterdata=spdata(spdata.treatmentCtrl_0_LN_1_HN_2 == treatmentcode(k),:);
            nexttile
            scatter(1-0.01*litterdata.M_M_0_, litterdata.ligninCInG./litterdata.ligninCInG(1), 50)
            hold on; plot([1,0],[0,1])
            title("i"+i+" j"+j+" k"+k)
        end
    end
end
%%
close all
par_=[];

for i= 1:2
    sitedata = data(data.SiteCode == i,:);
    spcode = unique(sitedata.SpeciesCode);
    for j =1:length(spcode)
        spdata=sitedata(sitedata.SpeciesCode == spcode(j),:);
        treatmentcode = unique(spdata.treatmentCtrl_0_LN_1_HN_2);
        for k =1:length(treatmentcode)
            litterdata=spdata(spdata.treatmentCtrl_0_LN_1_HN_2 == treatmentcode(k),:);
            C_2 = 0.01*litterdata.M_M_0_(2);
            if (C_2 < 0.7)
                tsrt = 2;
                temp = 1;
            else
                tsrt = 1;
                temp = 0;
            end

            AISC0 = litterdata.ligninCInG(tsrt)./litterdata.CInG(tsrt);
            CN0 = litterdata.CN0(tsrt);
            obs_data=[];
            obs_data.tobs  = (litterdata.month(tsrt:end) -litterdata.month(tsrt)) *30;
            obs_data.Ct_obs = litterdata.CInG(tsrt:end);
            obs_data.Co_obs  = aromatic_fraction_inAIS(litterdata.ligninCInG(tsrt:end));
            figure(1);clf
            scatter(obs_data.tobs,obs_data.Ct_obs./obs_data.Ct_obs(1));hold on
            scatter(obs_data.tobs, obs_data.Co_obs./obs_data.Co_obs(1));

            CT_N = obs_data.Ct_obs./obs_data.Ct_obs(1);
            iddiff = [false;diff(CT_N)./CT_N(1:end-1) >0.1];
            obs_data.tobs(iddiff)=[];obs_data.Ct_obs(iddiff)=[];obs_data.Co_obs(iddiff)=[];

            plot(obs_data.tobs,obs_data.Ct_obs./obs_data.Ct_obs(1),'o-');hold on
            plot(obs_data.tobs, obs_data.Co_obs./obs_data.Co_obs(1),'o-');
            ylim([0,inf])
            title("i"+i+" j"+j+" k"+k)

            tnorm = obs_data.tobs;
%             f = fit(obs_data.tobs ,obs_data.Ct_obs./obs_data.Ct_obs(1),'exp2');
            final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
%             fun=@(x) f(x)-final_C*fraction_of_final_C_at_Terminal_time;
%             terminalTime = fsolve(fun,2000);
            f = fit(obs_data.tobs ,obs_data.Ct_obs./obs_data.Ct_obs(1),'exp1');
            decayk=f.b ;
            terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/decayk;
            n=terminalTime./obs_data.tobs(end);
            tt=0:tnorm(end)*n;
            
            param.CO_0=obs_data.Co_obs(1);
            param.CT_0=obs_data.Ct_obs(1);
            param.emax = emax_fun(CN0);
            init_guess = [-decayk*100, 0.2,0.1]; % [vh_max, mo,ro];
            fig=figure(2);clf
            [par,sol,~,~] = find_parameter(obs_data,...
                param,init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
            [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
            param.rsquare = r2;
            param.rmse = rmse;
            makeplot_state_space(sol, '',param,obs_data,g,fig)

            exportgraphics(fig, fig_path_all+study+"_"+sitedata.Site(1)+"_"+ ...
                spdata.Var3(1) +"_k"+k+".png",'Resolution',100)
            par_ = [par_;j,k, par,AISC0,CN0,AISC0*CN0,rmse,r2];
            disp("current simulation is i"+i+" j"+j+" k"+k )
%             pause
        end
    end
end
%%
save_fname = excel_path+ "fitted_par_"+study+".txt";
save(save_fname,"par_",'-ascii','-double','-tabs')