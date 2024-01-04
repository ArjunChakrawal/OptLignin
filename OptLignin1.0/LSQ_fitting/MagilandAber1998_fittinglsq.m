clear all
clc
close all
%%

set_up
study="MagilandAber1998";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";

%%
data=readtable('..\data/Magil and Aber 1998.xlsx','Range','A8:Q242');
%%
close all
par_=[];
fig = figure;
fig.Position = [100, 100, 1200, 800];
tiledlayout('flow');
tau = [];
Ctau=[];
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

            CT_N = obs_data.Ct_obs./obs_data.Ct_obs(1);
            iddiff = [false;diff(CT_N)./CT_N(1:end-1) >0.1];
            obs_data.tobs(iddiff)=[];obs_data.Ct_obs(iddiff)=[];obs_data.Co_obs(iddiff)=[];


            param.emax = emax_fun(CN0);

            init_guess = [0.03	0.03 0.2	50];
            [par, sol] = find_parameter_lsq(obs_data, param, ...
                init_guess, @ysim_state_space_lsq)  ;
            [r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);

            nexttile
            fit_plotting(sol, obs_data);
            title(study+" "+sitedata.Site(1)+" "+ ...
                spdata.Var3(1) +" k"+k)
            CT=sol.y(1,:);Co=sol.y(2,:);

            dcodt = diff(sol.y(2,:))./diff(sol.x);
            neg_idx  = find(dcodt<0, 1);
            if(isempty(neg_idx))
                neg_idx=length(sol.x);
            end
            tau=sol.x(neg_idx);
            CTN=CT./CT(1); Ctau=CTN(neg_idx);


            par_ = [par_;par,AISC0,CN0,AISC0*CN0,...
                r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];
            disp("current simulation is i"+i+" j"+j+" k"+k )
        end
    end
end
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)





