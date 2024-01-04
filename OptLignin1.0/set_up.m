
param.numiter=50;
lw = 3;
aromatic_fraction_inAIS = @(AIS)0.2 * AIS;
szp = 50;
LC = linspecer(4);
scenario_name = "exp";

switch scenario_name
    case "exp"
        a=0.1537;
        b = 2;
        g = @(l, a, b) exp(-(l / a).^b);
end

fraction_of_final_C_at_Terminal_time = 0.5;
fraction_of_C_in_AUR = 0.6;
initial_fraction_of_C = 0.5;

vomax = 0.1;
emax = 0.3;
emax_fun = @(CN)min(6.25*CN.^-0.77, 0.4);
f = 1;
param.vomax = vomax;
param.f = f;
param.a = a;
param.b = b;
param.g=g;
param.mo = 0.1;
param.vo_thres = 0.05;
vo_thres=0.05;
