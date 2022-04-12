function [theta,mdlfit] = static_glm_adaptive(s_pre,s_post,fit,varargin)% 1- loading the spiking activity

dt   = 0.001; % binsize (don't change it)

% defaults
% delay=[100 100]/(dt*1000);
% nfilt = [5 5];
% NumLambda = 20;
% numCV = 8;
% toleranceValue=.1;

if (~isempty(varargin))
    c = 1 ;
    % user defined
    while c <= length(varargin)
        switch varargin{c}
            case {'delay'}
                delay = varargin{c+1};
            case {'nfilt'}
                nfilt = varargin{c+1};
            case {'NumLambda'}
                NumLambda = varargin{c+1};
            case {'numCV'}
                numCV = varargin{c+1};
            case {'toleranceValue'}
                toleranceValue = varargin{c+1};
        end % switch
        c = c + 2;
    end % for
end % if

% 2- create basis functions for coupling and post-spike history filter
% Fc = orth(getBasis('rcos',nfilt(1),delay(1),20,0)')';


Xc = getX(s_pre,fit.Fc,0,1,0)';
Xc = circshift(Xc,-floor(length(fit.Fc)/2));

Xs = getX(s_pre,fit.Fsmooth,0,1,0)';
Xs = circshift(Xs,-floor(size(fit.Fsmooth,2)/2));

x0 = linspace(0,1,1/dt);hist_tau = .01;
fit.Fh = [0 exp(-x0/hist_tau)];
Xh = getX(s_post,fit.Fh,0,1,0)';


% 6- glm-static
 % check if the T/dt is an integer
dev=0;
tic
fprintf('GLM-static...')
Hist = Xh*fit.hist_beta;
offset = Hist+log(dt);
bta = glmfit([Xc Xs],s_post,'poisson','offset',offset);
w = bta(2);
Coup = Xc*w;
Slow_eff = Xs*bta(3:end);
b0 = [ones(length(s_pre),1)*bta(1)];

% figure,
for i = 2%:40 % could change the maximum number of iteration
    %adaptive filter
    W0 = eye(1)*0.1;
    X = [ones(length(s_pre),1) ];
    a_offset = Hist+Coup+Slow_eff;
    [b,W_adapt,~] = adaptivefilter(s_post,X,b0(1,:)',W0,fit.F,fit.Q,a_offset,1);
%     plot(b);hold on;
    
    % stp model
    STP_offset = b'+ Hist + log(dt); % static part of the lambda
    clear stats
    [alph,dev(i),stats] = glmfit([Xc*w Xs],s_post','poisson','Offset',STP_offset);
    stats.Deviance = dev(i);
    
    % update parameters
    b0 = [b'+alph(1)];    
    w = w*alph(2);
    Coup = Xc*w;
    Slow_eff = Xs*alph(end-size(Xs,2)+1:end);
    
    fprintf('loop %d: Dev difference: %04.01f in %02.01f \n',i,(dev(i)-dev(i-1)),toc);
    if i>3 && dev(i) - dev(i-1)>0;break;end
    if abs(dev(i)-dev(i-1))<toleranceValue ;break;end % could change the convergence limit
     
end
yhat_static = exp(b0 + Coup + Slow_eff + Hist)*dt;
llhd_static_all = -yhat_static + log((yhat_static+(yhat_static==0))).*(s_post'); 

theta.w = w;
theta.b0 = b';
theta.alph = alph;
theta.bta = bta;
theta.stats = stats;
theta.fit = fit;

mdlfit.yhat = yhat_static;
mdlfit.yhat0 = yhat0_static;
mdlfit.llhd_all = sum(llhd_static_all);
mdlfit.Coup = Coup;
mdlfit.Slow_effect = Slow_eff;
mdlfit.Hist = Hist;

