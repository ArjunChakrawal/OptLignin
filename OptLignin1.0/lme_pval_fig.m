function lme_pval_fig(mdl,ax)
axes(ax)
p = mdl.Coefficients.pValue;
e = mdl.Coefficients.Estimate;

map = [212, 17, 89; ...
    mean([212, 17, 89; 255, 255, 255], 1); ...
    255, 255, 255; ...
    mean([26, 133, 255; 255, 255, 255], 1); ...
    26, 133, 255] / 255;



sz=size(p);

pval = 0.95; % significance threshold for all analyses
z = e./abs(e);
zz = p;
zz(p>1-pval) = 0;
zz(p>1-pval & p<0.1) = 0.5;
zz(p<=1-pval) = 1;
zzz = zz.*z;
image(zzz(2:length(p))','CDataMapping','scaled'); hold on
colormap(map)


for col=2:length(p)
    text(col-1, 1, "\beta="+num2str(e(col),2),'HorizontalAlignment','center');
    text(col-1, 0.75, "\itpval="+num2str(p(col),2),'HorizontalAlignment','center');
end

xticks(1:length(p)-1);
yticks([1]);
xticklabels(mdl.CoefficientNames(2:end))
yticklabels(mdl.ResponseName)
str1="\it AIC="+mdl.ModelCriterion.AIC;
str2="\it r^2="+mdl.Rsquared.Ordinary;
title(str1+" "+str2)