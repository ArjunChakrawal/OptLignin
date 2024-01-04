clear all
clc
close all
LC=linspecer(4);

%
% scenario_name="without_g_func";
% scenario_name="with_g_func";

scenario_name="exp_g";

lw = 3;
col = [[255, 69, 0] / 255; [30, 144, 255] / 255; [0, 0, 0]];
aromatic_fraction_inAIS = @(AIS)0.25*AIS;
A=0.5;B=0.05;
a = A*(0.25);b= B*(0.25);

L = 0:0.01:1;
l=aromatic_fraction_inAIS(L);
figure
plot(l, 1-1./(1 + exp(-(l - a)./b)), 'LineWidth', 2, 'DisplayName', "sigmoid"); hold on
plot(l, exp(-(l/0.1).^2),'LineWidth',2, 'DisplayName', "exponential")
ylabel("Factor reducing carbohydrate uptake")
xlabel('Lignocellulose index = Co/(Ch+Co)');
lh = legend("show");
lh.Box = 'off';
lh.FontSize = 12;
lh.Location = "best";
set(gca, 'FontSize', 13, 'Color', 'w')
% Read_data
Read_data
%
% close all
datacode = unique(berg_data.DatasetCode);
i = 1;
id = find(berg_data.DatasetCode == i);
dataset = berg_data(id(1):id(end), :);
spcode = unique(dataset.SpeciesCode);
j = 1;
idsp = find(dataset.SpeciesCode == spcode(j));
decdata2 = dataset(idsp(1):idsp(end), :);
id = decdata2.NMg_g < 0;
decdata2(id, :) = [];
total_C = (1-decdata2.MassLoss_*0.01)*initial_fraction_of_C; % gC/ bag
amount_AUR_C = (decdata2.AISMg_gInitialLitter*0.001*0.5).*(1-decdata2.MassLoss_(1)*0.01);  % gligC/ bag
aromaticC = aromatic_fraction_inAIS(amount_AUR_C); % true lignin


C_2 = decdata2.gC_gCInitial(2);
N_2 = decdata2.gN_gNInitial(2);
if (N_2 < 0.9 && N_2 < C_2 || C_2 < 0.7)
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

figure;
tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
nexttile
scatter(obs_data.tobs,obs_data.Ct_obs,'DisplayName','TC'); hold on
scatter(obs_data.tobs,obs_data.Co_obs,'DisplayName','AromaticC'); hold on
xlabel("time"); ylabel('gC/bag')

nexttile
scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1)); hold on
xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
% Fixed parameters
% close all
fig=figure;

param.emax = 0.3;

g = @(L,a,b)1 - 1./(1 + exp(-(L-a)./b));
% g = @(l,a,b)1 ;
% g = @(l, a, b) exp(-(l/0.1).^2);
init_guess = [0.0005 0.003	0.05	100]; % [vh_max, mo,ro]; 

[ysim,sol]= ysim_state_space_least_square(init_guess, param, obs_data,g,a,b);

figure
set(gcf, 'Color','w')
plot(sol.x./365,sol.y./sol.y(:,1)); hold on 
scatter(obs_data.tobs./365, obs_data.Ct_obs./obs_data.Ct_obs(1),...
    30,'DisplayName',"CTobs",'MarkerFaceColor',LC(1,:) ); 
try
    scatter(obs_data.tobs_Co./365, obs_data.Co_obs./obs_data.Co_obs(1), ...
        30, 'DisplayName',"COobs",'MarkerFaceColor',LC(2,:) ); hold on
catch
    scatter(obs_data.tobs./365, obs_data.Co_obs./obs_data.Co_obs(1), ...
        30, 'DisplayName',"COobs",'MarkerFaceColor',LC(2,:)); hold on
end

%
ydata = [obs_data.Ct_obs./obs_data.Ct_obs(1);...
obs_data.Co_obs./obs_data.Co_obs(1)];

fun = @(est_par, t) ysim_state_space_least_square(est_par, param, obs_data,g,a,b);

% ysim = fun(init_guess, obs_data.tobs);

options1 = optimoptions('lsqcurvefit', 'Display', 'final', ...
         'PlotFcn', {@optimplotx, @optimplotfirstorderopt, @optimplotresnorm}, ...
    'MaxFunEvals', 2000, 'MaxIter', 1000, 'TolFun', 1e-12, 'TolX', 1e-12,...
    'FiniteDifferenceStepSize', 1e-9, 'FiniteDifferenceType', 'forward');
lb = [0.00001,0.0001, 0.0001,1]; % [vomax,vh_max, mo,ro]
ub = [0.01,0.01, 0.5,500];

problem = createOptimProblem('lsqcurvefit',...
    'objective',fun,...
    'xdata',obs_data.tobs,'ydata',ydata,...
    'x0',init_guess,...
    'lb',lb,'ub',ub,...
    'options',options1);

[par, resnorm,~, ~, ~, ~, ~] = lsqcurvefit(problem);


% resn=[];parest=[];normres=[];
% fac=[1];
% for i =1:length(fac)
%     problem.x0(1)=init_guess(1)*fac(i);
%     [par, resnorm,~, ~, ~, ~, ~] = lsqcurvefit(problem);
%     [ysim, sol] = fun(par, obs_data.tobs);
%     rmse = sqrt(mean(ysim-ydata).^2);
%     resn=[resn;rmse];
%     parest=[parest;par];
%     normres=[normres;resnorm];
% end

[ysim,sol]= ysim_state_space_least_square(par, param, obs_data,g,a,b);

figure
set(gcf, 'Color','w')
plot(sol.x./365,sol.y./sol.y(:,1)); hold on 
scatter(obs_data.tobs./365, obs_data.Ct_obs./obs_data.Ct_obs(1),...
    30,'DisplayName',"CTobs",'MarkerFaceColor',LC(1,:) ); 
try
    scatter(obs_data.tobs_Co./365, obs_data.Co_obs./obs_data.Co_obs(1), ...
        30, 'DisplayName',"COobs",'MarkerFaceColor',LC(2,:) ); hold on
catch
    scatter(obs_data.tobs./365, obs_data.Co_obs./obs_data.Co_obs(1), ...
        30, 'DisplayName',"COobs",'MarkerFaceColor',LC(2,:)); hold on
end
figure;
plot(sol.y(1,:),sol.y(2,:))
%% Berg fitting
init_guess = [0.003, 0.1,90]; % [vh_max, mo,ro]; with parabolic g
ix=0;
fig=figure;
par_berg=[];
for i = 1:length(datacode)
    id = find(berg_data.DatasetCode == i);
    dataset = berg_data(id(1):id(end), :);
    spcode = unique(dataset.SpeciesCode);
    for j = 1:length(spcode)
        ix = ix + 1;
        idsp = find(dataset.SpeciesCode == spcode(j));
        decdata2 = dataset(idsp(1):idsp(end), :);
        id = decdata2.NMg_g < 0;
        decdata2(id, :) = [];
        total_C = (1-decdata2.MassLoss_*0.01)*initial_fraction_of_C; % gC/g litter
        amount_AUR_C = (decdata2.AISMg_gInitialLitter*0.001*0.5).*(1-decdata2.MassLoss_(1)*0.01);  % gligC/g litter
        aromaticC = aromatic_fraction_inAIS(amount_AUR_C); % true lignin

        C_2 = decdata2.gC_gCInitial(2);
        N_2 = decdata2.gN_gNInitial(2);
        if (N_2 < 0.9 && N_2 < C_2 || C_2 < 0.7)
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
        % obs_data.Co_obs = amount_AIS_C(tsrt:end);
%         obs_data.T = decdata2.Time_days_(end);
        param.CO_0=obs_data.Co_obs(1);
        param.CT_0=obs_data.Ct_obs(1);

        [par,sol,rmse,rsquare] = find_parameter( obs_data,...
            param,init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);

        AISC0=decdata2.gLigninC_gC(tsrt);
        CN0=decdata2.gCInitial_gNInitial(tsrt);
        LN0=decdata2.gligC0_gN0(tsrt);
        par_berg = [par_berg;i,j, par,AISC0,CN0,LN0,rmse,rsquare];
        i
        j
    end
end
save_fname = scenario_name+ "_Berg.txt";
save(save_fname,"par_berg",'-ascii','-double','-tabs')

%% Hirobe fitting
close all
Initial_litter_mass = 5; % grams
datacode = unique(organic_hirobe.Species);
init_guess = [0.03, 0.15,10]; % [vh_max, mo,ro]; with parabolic g for berg
fig=figure;

par_Hirobe=[];
for i = 1:length(datacode)
    id = find(Hirobe04_data.Species == datacode(i));
    ddata = Hirobe04_data(id, :);
    obs_data=[];
    obs_data.tobs  = ddata.Time*30;
    Nt_obs = Initial_litter_mass.*(0.01 .* ddata.weightRemain_ )...
        .* (ddata.N_ .* 0.01);
    Ct_obs = Initial_litter_mass.*(0.01 .* ddata.weightRemain_ )...
        .* (ddata.C_ .* 0.01); % g litter/ bag * fraction of litter*fraction of C = (gC/g bag)
    obs_data.Ct_obs = Ct_obs;

    datacode = unique(organic_hirobe.Species);
    idL = find(organic_hirobe.Species == datacode(i));
    data_Lig_Carb = organic_hirobe(idL, :);
    obs_data.tobs_Co = data_Lig_Carb.Period.*30;

    for j = 1:length(data_Lig_Carb.Period)
        idtemp(j) = find(ddata.Time == data_Lig_Carb.Period(j));
    end
    % g litter/ bag * fraction of litter*fraction of lignin C (g lignin C/ g litter)= g lignin C/ bag
    AIS_C = Initial_litter_mass.*0.01 * ddata.weightRemain_(idtemp) .* ...
        data_Lig_Carb.Lignin_mg_g_ .* 0.001*fraction_of_C_in_AUR;
    aromaticC = aromatic_fraction_inAIS(AIS_C); % true lignin
    obs_data.Co_obs=aromaticC;
    param.CO_0=obs_data.Co_obs(1);
    param.CT_0=obs_data.Ct_obs(1);

    [par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
        init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);

    par_Hirobe = [par_Hirobe;i, par,AIS_C(1)/Ct_obs(1),Ct_obs(1)/Nt_obs(1),...
        AIS_C(1)/Nt_obs(1),rmse,rsquare];

end
save_fname = scenario_name+ "_par_Hirobe.txt";
save(save_fname,"par_Hirobe",'-ascii','-double','-tabs')

%% Ososno et al. 2004
close all

init_guess = [0.03, 0.15,100]; % [vh_max, mo,ro]; with parabolic g for berg
sz =size(amount_C);
par_Osono=[];
fig=figure;
for i = 1:sz(2)
    gC_gCInitial = amount_C(:, i)./amount_C(1, i);
    C_2 = gC_gCInitial(2);
    if ( C_2 < 0.7)
        tsrt = 2;
        temp = 1;
    else
        tsrt = 1;
        temp = 0;
    end
    tt=osono_t(tsrt:end)-osono_t(tsrt);
    obs_data = [];
    obs_data.tobs = tt;
    obs_data.Ct_obs = amount_C(tsrt:end, i); %Remaining mass of litter (gC / bag)
    obs_data.Co_obs = aromatic_fraction_inAIS(amount_AUR_C_osono04(tsrt:end, i));
    param.CO_0=obs_data.Co_obs(1);
    param.CT_0=obs_data.Ct_obs(1);
    [par,sol,rmse,rsquare]  = find_parameter( obs_data, param,...
        par,g,fig,@makeplot_state_space,@ysim_state_space,n);

    par_Osono = [par_Osono;i, par,obs_data.Co_obs(1)/obs_data.Ct_obs(1) ...
        ,CN0_osono(i),LN0_osono(i),rmse,rsquare];
end
save_fname = scenario_name+ "_par_Osono.txt";
save(save_fname,"par_Osono",'-ascii','-double','-tabs')

%% Osono2017
close all
fig=figure;
init_guess = [0.03, 0.15,100]; % [vh_max, mo,ro]; with parabolic g for berg
par_Osono17=[];
numSpecies = 12;
for i = 1:numSpecies
    if (i == 12)
        ddata = Osono2017(iddata(i):length(Osono2017.Collection), :);
    else
        ddata = Osono2017(iddata(i):iddata(i+1)-1, :);
    end
    id_t = ddata.Collection;
    tt = [0, month(id_t(2:end))] .* 30;
    Time_days_ = tt';
    amount_C = ddata.MassG .* ddata.C_ .* 0.01; % g C/ bag
    amount_N = ddata.MassG .* ddata.N_ .* 0.01; % g N/ bag
    amount_AUR_C = fraction_of_C_in_AUR.*ddata.MassG .* ddata.AUR_ .* 0.01; % g lignin / bag

    gC_gCInitial = amount_C./amount_C(1);
    C_2 = gC_gCInitial(2);
    if ( C_2 < 0.7)
        tsrt = 2;
        temp = 1;
    else
        tsrt = 1;
        temp = 0;
    end
    tt=Time_days_(tsrt:end)-Time_days_(tsrt);
    obs_data = [];
    obs_data.tobs = tt;
    obs_data.Ct_obs = amount_C(tsrt:end); %Remaining mass of litter (gC / bag)
    obs_data.Co_obs = aromatic_fraction_inAIS(amount_AUR_C(tsrt:end));
    param.CO_0=obs_data.Co_obs(1);
    param.CT_0=obs_data.Ct_obs(1);

    [par,sol,rmse,rsquare]  = find_parameter( obs_data, param,...
        init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
    AISC0=amount_AUR_C(tsrt)/amount_C(tsrt);
    CN0=amount_C(tsrt)/amount_N(tsrt);
    LN0=amount_AUR_C(tsrt)/amount_N(tsrt);
    par_Osono17 = [par_Osono17;i, par,AISC0,CN0,LN0,rmse,rsquare];
end
save_fname = scenario_name+ "_par_Osono17.txt";
save(save_fname,"par_Osono17",'-ascii','-double','-tabs')




