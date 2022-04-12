function [theta,mdlfit] = stp_gblm_adaptive_v5(s_pre,s_post,fit,varargin)% 1- loading the spiking activity

dt   = 0.001; % binsize (don't change it)
L = length(s_pre);
isi = diff(find(s_pre>0)*dt); % interspike intervals

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

% 3- spline expansion
% Nq=8; % 8-20 suggested but the more baseline = overfitting
Nq = fit.Nq;
Bm = getBasis('rcos',Nq,1000,20,0);
Bm = padarray(Bm,[0 5000]);
Bm=Bm(:,5001:end);
isi_modi = isi;
isi_modi(isi>5) = 5;
Bm_dt = Bm(:,round(isi_modi*1000));
s=zeros(Nq,L);
for i =1: Nq
    s(i,s_pre>0)=[0 Bm_dt(i,:)];
end

% 4- transient effect of stp carried by exp function
% from splined spikes to exponential ( from "s" [Nq x T] --> "e" [Nq x T] )
% tau= 0.2; 
tau = fit.tau;% timeconstant for decay
x0 = linspace(0,1,1000);
kern_stp = exp(-x0/tau);
e = filter(kern_stp,1,s');
% e = circshift(e,7);


% 5- gblm
 % check if the T/dt is an integer
dev=0;
tic
Hist = Xh*fit.hist_beta;
offset = Hist+log(dt);
bta = glmfit([Xc Xs],s_post,'poisson','offset',offset);
w = bta(2);
wt = ones(L,1);
Coup_modi = Xc*w.*wt;
Slow_eff = Xs*bta(3:end);
b0 = [ones(length(s_pre),1)*bta(1)];

% figure,
for i = 2:40 % could change the maximum number of iteration
    %adaptive filter
    W0 = eye(1)*0.1;
    X = [ones(length(s_pre),1) ];
    a_offset = Hist+Coup_modi+Slow_eff;
    [b,~,~] = adaptivefilter(s_post,X,b0(1,:)',W0,fit.F,fit.Q,a_offset,1);
%     plot(b);hold on;
    
    % stp model
    STP_offset = b'+ Hist + log(dt); 
    W_short = repmat(Xc*w,1,Nq).*e; % plastic part of lambda  
    clear stats
    [alph,dev(i),stats] = glmfit([Xc*w W_short Xs],s_post','poisson','Offset',STP_offset);
    stats.Deviance = dev(i);
    
    % update parameters
    b0 = [b'+alph(1)];    
    w = w*alph(2);
    wt = 1 + e*(alph(3:end-size(Xs,2))/alph(2));
    Coup = Xc*w;
    Coup_modi = Xc*w.*wt;% modified coupling term
    Slow_eff = Xs*alph(end-size(Xs,2)+1:end);
    
    fprintf('loop %d: Dev difference: %04.01f in %02.01f \n',i,(dev(i)-dev(i-1)),toc);
    if i>3 && dev(i) - dev(i-1)>0;break;end
    if abs(dev(i)-dev(i-1))<toleranceValue ;break;end % could change the convergence limit
     
end
yhat = exp(b0 + Coup_modi + Slow_eff + Hist)*dt;
llhd_all = -yhat + log((yhat+(yhat==0))).*(s_post');  

% var-cov matrix
offset = b0 + Hist + Slow_eff + log(dt);
W_short = repmat(Xc*w,1,Nq).*e;
[~,~,stats2] = glmfit(W_short,s_post','poisson','Offset',offset,'constant','off');
se_modif_fxn = sqrt(diag(Bm'*stats2.covb*Bm));

% parameters & hyperparameters
theta.w = w;
theta.b0 = b';
theta.alph = alph;
theta.bta = bta;
theta.modif_fxn = alph(3:end-size(Xs,2))'/alph(2)*Bm;
theta.se_modif_fxn = se_modif_fxn';
theta.Bm = Bm;
theta.stats = stats;
theta.fit = fit;

%model fits
mdlfit.yhat = yhat;
mdlfit.yhat0 = yhat0;
mdlfit.llhd_all = sum(llhd_all);
mdlfit.Coup = Coup;
mdlfit.Coup_modi = Coup_modi;
mdlfit.Slow_effect = Slow_eff;
mdlfit.Hist = Hist;


