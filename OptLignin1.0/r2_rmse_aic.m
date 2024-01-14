function [r2, rmse,AIC, BIC] = r2_rmse_aic(ydata, ysim,p)

BICfun = @(n, p, ssr) n * log(ssr/n) + p*log(n);

res=ysim-ydata;
dev= mean(ydata) -ydata;
SST = sum(dev.^2);
SSR = sum(res.^2);
r2 = 1 - SSR / SST;
rmse=sqrt(mean(res.^2));

% now we are calculating BIC because AIC can inaccurate values if n - p - 1
% is zero
AIC =BICfun(length(ydata),p,SSR); 
% BIC = BICfun(length(ydata),p,SSR);
% AICfun = @(n, p, ssr) n * log(ssr/n) + 2 * n * p / (n - p - 1);
% AICfun = @(n, p, ssr) n * log(ssr/n) + 2  * p;

    function out=AICfun(n, p, ssr)
        if(n - p - 1~=0)
            out=n * log(ssr/n) + 2 * n * p / (n - p - 1);
        else
            out=n * log(ssr/n) + 2  * p;
        end
    end

end
