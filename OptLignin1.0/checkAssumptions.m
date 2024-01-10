function [lme,tbl] = checkAssumptions(data, target, pred_formula, fname, varargin)
% Fit and performs diagnostic checks for a linear regression model

p = [];
e = [];
r2=[];AIC=[];SEe=[];
lme=cell(1,length(target));
tbl=cell(1,length(target));

for i = 1:length(target)
    formula =  [target{i}, pred_formula];
    switch fname
        case 'lme'
            mdl = fitlme(data, formula, varargin{:});
        case 'lm'
            mdl = fitlm(data,formula, varargin{:});
    end
    % Check assumptions

    % 1. Linearity check
    fig = figure;
    fig.Position(3:4) = [ 850, 600];
    fig.Color = 'w';
    tiledlayout('flow')
    try
        plot(mdl);
        nexttile
        plot(mdl);
        title('Linearity check');
        legend("r^2 = " +mdl.Rsquared.Ordinary)
    catch
        nexttile
        scatter(data.(mdl.ResponseName), fitted(mdl))
        xlabel('observation'); ylabel('fitted values')
        title('Linearity check');
        legend("r^2 = " +mdl.Rsquared.Ordinary)
    end
    % 2. Normality check
    nexttile
    % Extract residuals
    residuals = mdl.Residuals.Raw;
    qqplot(residuals);
    title("Normality Check: residual QQ plot of model")

%     nexttile
%     plotResiduals(mdl, 'probability');
    % 3. Constant variance check
    nexttile
    plotResiduals(mdl, 'fitted');
    title('Constant variance check');

    % 4. Independence of errors check
    nexttile
    plotResiduals(mdl, 'lagged');
    title('Independence of errors check');
    % Diagnostic Plots
    try
        nexttile
        plotDiagnostics(mdl,'cookd')     % Look for points with large Cook’s distance.
    catch
    end
    B = randomEffects(mdl);
    nexttile
    qqplot(B);
    title("Normality of random effects")

    str="model:" +mdl.Formula.char;
    sgtitle(str,'fontsize',11)

    p = [p, mdl.Coefficients.pValue];
    e = [e, mdl.Coefficients.Estimate];
    SEe = [SEe, mdl.Coefficients.SE];
    r2= [r2,mdl.Rsquared.Ordinary];
    AIC= [AIC,mdl.ModelCriterion.AIC];
    lme{i}=mdl;
    tbl{i}= anova(mdl);
%     exportgraphics(fig, "results\Diagnostic_plot_"+target{i}+".png", Resolution=300)
end
e;
p;


%% plot R2 and AIC
fig=figure;tiledlayout('flow');nexttile
fig.Position(3:4)=[ 600   200];
set(gcf, 'Color', 'w')

try
    subplot(121)
    bar(categorical(target),r2);
    labels1 = string(r2);
    text(categorical(target), ...
        r2,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    ylabel('\it R^2')
    ylim([0,1])
    subplot(122)
    bar(categorical(target),AIC);
    labels1 = string(AIC,'%.3f');
    text(categorical(target), ...
        AIC,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    ylabel('\it AIC')
    %     ylim([0,inf])
catch
end
% sgtitle("Y"+ pred_formula,'fontsize',11)
% exportgraphics(fig, 'results\R2_AIC_plot.png', Resolution=300)



%% image plot for estimates and their significance
rsquared_values = readtable('rsquared_values.csv');
r2_marg = rsquared_values.marginal_MarginalR2;
%%

fig=figure;tiledlayout('flow')
fig.Position(3:4)=[1000, 307];
set(gcf, 'Color', 'w')
nexttile
pval = 0.95; % significance threshold for all analyses
est = e ./ abs(e); % vector of postive:+1 or negative:-1 effect
% zz map pvalue significant:1, marginally significant:0.5, and not significant:0
zz = p;
zz(p <= 1-pval) = 1;
zz(p > 1-pval) = 0;
zz(p > 1-pval & p < 0.1) = 0.5;

zzz = zz .* est;
sz_zzz=size(zzz);
im=zzz(2:sz_zzz(1), :)';
sz_im=size(im);
im=[im,zeros(sz_im(1),1)];
image(im, 'CDataMapping', 'scaled'); hold on
% heatmap(zzz(2:length(p), :)')
% Define the color map
cmap = [255, 102, 102; % red
    255, 255, 255; % white
    0, 153, 153]/255; % blue

% Create a colormap with 256 levels
numLevels = 5;
customCmap = interp1(linspace(-1, 1, 3), cmap, linspace(-1, 1, numLevels));
colormap(customCmap)

% e= [e;r2_marg'];
sz=size(p);
for j = 1:length(target)
    for i = 2:sz(1)
        %         if zzz(i, j) ~= 0
        if  p(i,j)>0 && p(i,j)<0.001
            str=[num2str(e(i, j), 2),'±',num2str(SEe(i, j), 2)]+" ***";
        elseif p(i,j)>0.001 && p(i,j)<0.01
            str=[num2str(e(i, j), 2),'±',num2str(SEe(i, j), 2)]+" **";
        elseif p(i,j)>0.01 && p(i,j)<0.05
            str=[num2str(e(i, j), 2),'±',num2str(SEe(i, j), 2)]+" *";
        elseif p(i,j)>0.05 && p(i,j)<0.1
            str=[num2str(e(i, j), 2),'±',num2str(SEe(i, j), 2)]+" .";
        else
            str=[num2str(e(i, j), 2),'±',num2str(SEe(i, j), 2)]+" n.s.";
        end
        text(i-1, j, str, 'HorizontalAlignment', 'center',...
            'FontWeight','normal',FontSize=11,FontName='Serif',Color=[1 1 1]*0.1);        
    end
end

for i = 1:length(target)
    str = ['(',num2str(rsquared_values.marginal_MarginalR2(i), 2),...
    ' , ',num2str(rsquared_values.conditional_ConditionalR2(i), 2),')'];
    text(7,i, str, 'HorizontalAlignment', 'center',...
            'FontWeight','normal',FontSize=11,FontName='Serif',Color=[1 1 1]*0.1); 
end

sz_im=size(im);
xticks(1:sz_im(2));
yticks(1:length(target));
xticklabels([mdl.CoefficientNames(2:end),"({\itr^2_{marg} , r^2_{cond}})"])
yticklabels(target)
title("Y"+ pred_formula,'fontsize',11)
% str='Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1';
% p=scatter(nan,nan);
% lh=legend(p,str,'Location','southoutside');
% lh.Box="off";
ax=gca;
ax.TickLength=[0,0];
if(length(target)>1)
    h = pixelgrid;
    h.Children(1).Color = [1, 1, 1] * 0.8;
    h.Children(1).LineWidth = 0.5;
    h.Children(2).Color = [1, 1, 1] * 0.8;
    h.Children(2).LineWidth = 0.5;
end



end
