clear all
clc
close all
%%

set_up
study ="Tu2014";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
%%

data=readtable('..\data/Tu_et_al2014.xlsx','Range','B16:O360');
CN0=table2array(readtable('..\data/Tu_et_al2014.xlsx','Range','K2:K11'));
AISC0=table2array(readtable('..\data/Tu_et_al2014.xlsx','Range','N2:N11'));
SP=readtable('..\data/Tu_et_al2014.xlsx','Range','A1:A11');
figure;tiledlayout("flow")
for i =1:10
    nexttile
    tempdata1=data(data.SpeciesCode == i,:);
    gscatter(tempdata1.C_C_0_, tempdata1.lignin_lignin_0_, tempdata1.treatment)
    hold on; plot(0:1,0:1)
    title(SP.LitterSubstrate{i})
end
%

treatementName = unique(data.treatment);
SPname = unique(data.SpeciesCode);
%%
par_=[];
fig = figure;
fig.Position = [100, 100, 1200, 800];
tiledlayout('flow');
tau = [];
Ctau=[];

for i =1:length(SPname)
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

        CT_N = obs_data.Ct_obs./obs_data.Ct_obs(1);
        iddiff = [false;diff(CT_N)./CT_N(1:end-1) >0.1];
        obs_data.tobs(iddiff)=[];obs_data.Ct_obs(iddiff)=[];obs_data.Co_obs(iddiff)=[];

        param.emax = emax_fun(CN0(i));
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

        par_ = [par_;par,AISC0(i),CN0(i),AISC0(i)*CN0(i),...
            r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];
        disp("current simulation is i ="+ i + " j="+ j )

    end

end
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)





