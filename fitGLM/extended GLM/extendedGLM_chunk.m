function model_fits = extendedGLM_chunk(CCG, X, hyperparameter,Info,model_fits_all)

% This function fits all the cross-correlogrms from one presynaptic neuron.
%
% model_fits is a structure with 4 fields:
%
% -glm: the model fitting results of the slow model
% ---yhat: fitted CCG
% ---llh: log likelihood
% ---b: parameter
% -exc: the model fitting results of the excitatory full model (stage 2,
% with latency constraints)
% ---eta: [eta_w eta_dt eta_tau], eta_dt is calculated based on the
% estimation
% ---v: estimated "conduction velocity", dt = distance/v(1)+v(2)
%
% -inh: the model fitting results of the inhibitory full model (stage 2,
% with latency constraints). The fields have the same meaning with field exc.
%
% -stage1: the model fitting results of the full model (stage 1), with a
% field for exc and inh stage 1 respectively.



%%
interval = Info.interval;
binsize = Info.binsize;
tau0 = hyperparameter.tau0;
eta_w = hyperparameter.eta_w;
eta_tau = hyperparameter.eta_tau;

NN = length(CCG);
t = linspace(-interval/2,interval/2,interval/binsize+2);
t = t+mean(diff(t))/2;
t = t(1:interval/binsize+1);

XX = [X;ones(1,interval/binsize+1)];

distance = zeros(NN,1);


eta = [eta_w 0 eta_tau];


%% excitatory model
% fprintf('exc model...')
% slow model & full model stage 1

% v0 = [0;0];
% b_glm = nan(NN,7);
% yhat_glm = nan(NN,length(t));
yhat_slow = nan(NN,length(t));
llh_glm = nan(NN,1);

b_s1 = nan(NN,2);
covb = cell(NN,1);
yhat_s1 = nan(NN,length(t));
llh_s1 = nan(NN,1);



syn =  model_fits_all.exc.syn;
basis = model_fits_all.exc.b(1:size(XX,1))*XX;


for post= 1:NN
    
    if isempty(CCG{post})
        continue
    end
    
    y = CCG{post};
    
    % full model Stage 1 (without latency constraints):
    
    [b,~,stats] = glmfit([basis',syn(1:length(y))'],y','poisson','constant','off');
    covb{post} = stats.covb;
    b_s1(post,:) = b';
    yhat_s1(post,:) = exp(b'*[basis',syn(1:length(y))']');
    yhat_slow(post,:) = exp(b(1:size(basis,1))'*basis);
    lam = yhat_s1(post,:);
    llh_s1(post) = nansum((y.*log(lam+(lam==0))-lam));
    
end

llr = llh_s1 - llh_glm;


%% output

% model_fits.glm.yhat = yhat_glm;
% model_fits.glm.llh = llh_glm;
% model_fits.glm.b = b_glm;

model_fits.exc.yhat = yhat_s1;
model_fits.exc.yhat_slow = yhat_slow;
model_fits.exc.llh = llh_s1;
model_fits.exc.b = b_s1;
model_fits.exc.covb = covb;
model_fits.exc.eta = eta;
model_fits.exc.basis = [basis',syn(1:length(y))'];
