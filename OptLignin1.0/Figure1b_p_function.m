clear all;
close all;
clc

CN = readtable('data\Bonanomi 2011 NPH\C and N.xlsx', 'Sheet', 'CN', 'Range', "A1:F65");
NMR = readtable('data\Bonanomi 2011 NPH\C and N.xlsx', 'Sheet', 'NMR', 'Range', "A2:L66");

SP = unique(NMR.Species);
kCarb = [];
aromatic = [];

for i = 1:length(SP)
    temp1 = NMR(strcmp(SP(i), NMR.Species), :);
    temp2 = CN(strcmp(SP(i), NMR.Species), :);

    OandHaromatic = (temp1.O_substitutedAromaticC_phenolicAndO_arylC_ ...
        +temp1.H_AndC_substitutedAromaticC).* 0.01; % fraction of aroamtic in total C
    non_aromatic = (1-OandHaromatic) .* temp2.massReamainng_; % mass of nonaroamtic at each time point
    k30 = -(1 / 30) * log(non_aromatic(2)/non_aromatic(1)) * 365;
    kCarb = [kCarb; k30];
    aromatic = [aromatic; OandHaromatic(1)];

end

final_data = [ aromatic, kCarb ];


N = length(aromatic);
l = 0:0.001:0.5;
y_norm = kCarb ./ max(kCarb);

modelfun = @(b, x) b(2) .* exp(-(x / b(1)).^2);
opts = statset('nlinfit');

[bestbeta, R, J] = nlinfit(aromatic, kCarb, modelfun, [0.15, 14], opts);
ci = nlparci(bestbeta,R,'Jacobian',J);

[ypred, delta] = nlpredci(modelfun, l, bestbeta, R, 'Jacobian', J);
lower = ypred - delta;
upper = ypred + delta;

LC =linspecer(1);
figure;
yp = modelfun(bestbeta, aromatic) ./ bestbeta(2);
scatter(aromatic, y_norm,'ok','LineWidth',3); hold on
hold on
plot(l, ypred./bestbeta(2), 'Color',LC, 'LineWidth', 4)
plot(l, lower./bestbeta(2), '--','Color',LC, 'LineWidth', 2)
plot(l, upper./bestbeta(2), '--','Color',LC, 'LineWidth', 2)
% ylabel("$p(l)$",'interpreter', 'latex')
% p=plot(nan,nan, 'Color',LC, 'LineWidth', 2);
% lh=legend(p,"$ p=exp{\left(-\left(\frac{l}{a}\right)^2\right)}$",'Fontsize',16, ...
%     'interpreter', 'latex');
% lh.Box='off';
% xlabel('$\it l= C_O/(C_O+C_H)$ [Aromatic C \ total C$^{-1}$]','interpreter', 'latex');
set(gca, 'Fontsize',16, 'LineWidth', 0.45, ...
    'Xcolor', [1, 1, 1]*0, 'Ycolor', [1, 1, 1]*0)
set(gcf,'color','w')
xlim([0, 0.5]); ylim([0, 1])
rss = sum((y_norm - yp).^2);
sst = sum((y_norm - mean(y_norm)).^2);
r2 =  1 - rss / sst;
exportgraphics(gcf, 'results\Figure1b_p_function.png','Resolution',300)
