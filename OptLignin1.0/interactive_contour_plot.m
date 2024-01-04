
logAN0 = linspace(-1.5, 2, 50);
MAT = linspace(-2, 3, 50);
[logAN0, MAT] = meshgrid(logAN0, MAT);
logARC0 = logAN0 * 0 + -0.5; % LOW logARC0
logAN0_T = MAT .* logAN0;

f=figure;
slider = uicontrol('Style', 'slider', 'Position', [10 10 300 20], ...
    'Min', -2, 'Max', -0.5, 'Value', -1, ...
    'Callback', @(source,~) updateLogARC0(source, logAN0, MAT, tau_lme.Coefficients));
function updateLogARC0(source, logAN0, MAT, coefficients)
    logARC0Value = get(source, 'Value');
    logARC0 = logAN0 * 0 + logARC0Value;
    logAN0_T = MAT .* logAN0;
    Z = coefficients.Estimate(1) + ...
        coefficients.Estimate(2) .* MAT + ...
        coefficients.Estimate(3) .* logARC0 + ...
        coefficients.Estimate(4) .* logAN0 + ...
        coefficients.Estimate(5) .* logAN0_T;
    contourf(logAN0, MAT, Z, 'ShowText', 'on')
    caxis([-5,5])
    xlabel("log(AN0)");
    ylabel("MAT_S");
    title("Low logARC0 = "+logARC0Value);
    c = colorbar;
    hc = get(c, 'Title');
    set(hc, 'String', "\it \tau");
    colormap jet
end

