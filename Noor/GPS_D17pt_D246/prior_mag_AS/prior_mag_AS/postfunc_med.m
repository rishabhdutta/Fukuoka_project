function[posterior,logposterior,fracafter]=postfunc_med(pars,frac,varargin)
%logposterior needed to feed for the gibbssample program which 
%required to generate candidates in MH ('proprnd')

%pars=pars'; %changed for MH

load Hypo_FukUTM.mat

moofpara=size(pars,1);
params=pars';

LB = [   590   3725    605     3715    0   0   -150    -4.0     0     0   ]';
UB = [   610   3745    630     3735    30  20   -40     -0.5     0     0   ]';

logposterior=zeros(size(pars,1),1);  %zeros(1,size(pars,2));
posterior=zeros(size(pars,1),1);

for i=1:size(pars,1)        %size(pars,2)
    model=pars(i,:);
    model=model';
    [res]=OneFault_obj_end(model,varargin{:});
    
    leng = sqrt((model(3)-model(1)).^2 + (model(4)-model(2)).^2);
    wid = model(5);
    slip = sqrt(model(8).^2+model(9).^2);
    moment=2.5*leng.*wid.*abs(slip).*1e23;
    Mw=2/3*log10(moment)-10.7;
   
    t= (LB <= model & UB >= model);
    
    if(~isempty(find(t == 0,1)) )
        logposterior(i) = 1e-50;
        posterior(i)=0;
        
    else
        
        
        if (Mw>=6 && Mw<=7)
            
            %put a prior for moment magnitude
            normal=@(x,mu,sigmasquare)((1/sqrt(2*pi*sigmasquare))*exp(-0.5*((x-mu)/sqrt(sigmasquare)).^2));
            fraction=normal(Mw,6.6,frac);
            fraction=fraction/normal(6.6,6.6,frac);
             
            nlsmodel = convert_endnls(model');
            trian = Planar_to_tri(nlsmodel,0);
            after = prior_aftershocks(trian,Latutm,Lonutm,depth);
            fracafter = exp(-after^2*1.5625e-4);   %fracafter = exp(-after^2*.004);
            
            if frac > 100
                fraction =1;
            end
            logposterior(i)=-0.5*res;  %logposterior(i)=-0.5*res
            posterior(i)=exp(logposterior(i))*fraction*fracafter;
            
            
        else
            logposterior(i) = 1e-50;
            posterior(i)=0;
        end
    end
        
end

end
