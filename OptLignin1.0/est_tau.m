function [tau_L,C_tau,avgvo,init_vo,max_vo,vo_at_tau,obj,sol] =est_tau(param,ARC0_loop,ro_loop,...
    vh_loop,pnamef,vh,ro,ARC0,mo,tstep)

% a=param.a;
% b=param.b;
g=param.g;
vo_thres=param.vo_thres;

param.CO_0 = param.CT_0 .*ARC0_loop;
[ocp, ~] = opt_con(param, g, [vh_loop, param.mo, ro_loop], param.TT);
ocp.solve('controlIntervals', 100, ...
    'collocationPoints', 'legendre', ... %
    'polynomialDegree', 3, ...
    'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
    );
switch pnamef
    case 'vh'
        nloop=length(vh);
    case 'ro'
        nloop=length(ro);
    case 'ARC0'
        nloop=length(ARC0);
    case 'mo'
        nloop=length(mo);
end
tau_L = zeros(nloop, 1);
C_tau=tau_L;
avgvo = tau_L;
init_vo= avgvo;
max_vo= avgvo;vo_at_tau=avgvo;
obj= zeros(nloop, 1);


for i = 1:nloop
    switch pnamef
        case 'vh'
            if(tstep(i)~=100)
                ocp.solve('controlIntervals', tstep(i), ...
                    'collocationPoints', 'legendre', ... %
                    'polynomialDegree', 3, ...
                    'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
                    );
            end
            par=[vh(i), param.mo, ro_loop];
            IC = [param.CT_0, param.CT_0 .* ARC0_loop];
            ocp.Parameter(1:3).setBound(par)
            ocp.State(1:2).setInitial(IC)
            ocp.parameterize();
            sol = ocp.optimize( ...
                'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
                );
        case 'ro'
            par=[vh_loop, param.mo, ro(i)];
            IC = [param.CT_0, param.CT_0 .* ARC0_loop];
            ocp.Parameter(1:3).setBound(par)
            ocp.State(1:2).setInitial(IC)
            ocp.parameterize();
            sol = ocp.optimize( ...
                'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
                );

        case 'ARC0'
            par=[vh_loop, param.mo, ro_loop];
            IC = [param.CT_0, param.CT_0 .* ARC0(i)];
            ocp.Parameter(1:3).setBound(par)
            ocp.State(1:2).setInitial(IC)
            ocp.parameterize();
            sol = ocp.optimize( ...
                'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
                );
        case 'mo'
            par=[vh_loop, mo(i), ro_loop];
            IC = [param.CT_0, param.CT_0 .* ARC0_loop];
            ocp.Parameter(1:3).setBound(par)
            ocp.State(1:2).setInitial(IC)
            ocp.parameterize();
            sol = ocp.optimize( ...
                'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
                );
    end

    %     emax=param.emax;
    time = sol.NumericalResults.Independent;
    vo = sol.NumericalResults.Control(1, :);
    CT = sol.NumericalResults.State(1, :);
    Co = sol.NumericalResults.State(2, :);
    %     vh_sol = sol.NumericalResults.Parameter(1);
    %     ro_sol = sol.NumericalResults.Parameter(3);
    %     f=param.f;
    %     L = Co ./ CT;
    %     gL = g(L, a, b);
    %     Dh = vh_sol .* gL .* (CT - Co);
    %     Do = vo .* Co;
    %     e = emax - ro_sol * vo;
    %     G = e .* (Dh + f .* Do);

    obj(i)=abs(sol.NumericalResults.Objective);
    %     if(obj(i)-abs(sol.NumericalResults.Objective)>1e-4)
    %         error("High numerical error in computation")
    %     end
    dcodt = diff(Co)./diff(time);
    neg_idx  = find(dcodt<0, 1);

    voN = vo ./ max(vo);
    if(isempty(neg_idx ))
        id=length(time);
    else
        id = find(voN > vo_thres);
    end
    tau_L(i) = time(id(1)) ;
    C_tau(i) = CT(id(1))/CT(1);
    avgvo(i) = mean(vo);
    init_vo(i) = vo(1);
    max_vo(i)=max(vo);
    vo_at_tau(i)=vo(id(1));

end
end