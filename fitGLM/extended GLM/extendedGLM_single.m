function model_fits = extendedGLM_single(CCG, ACG, X, hyperparameter,Info)

% This function fits all the cross-correlogrms from one presynaptic neuron.
% input:
% -CCG: cross-correlations
% -ACG: auto-correlation for the presynaptic neuron
% -X: basis functions
% -hyperparameter: the hyperparameters for w and tau
% -Info: the interval and binsize of the CCG
%
% output:
% model_fits is a structure with 3 fields:
% -glm: the model fitting results of the slow model
% ---yhat: fitted CCG
% ---llh: log likelihood
% ---b: parameter
% -exc: the model fitting results of the excitatory full model 
% ---yhat: fitted CCG
% ---yhat_slow: fitted CCG with w set to 0
% ---llh: log likelihood
% ---b: parameter
% -inh: the model fitting results of the inhibitory full model. The fields have the same meaning with field exc.
%



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
fprintf('exc model...')
% slow model & full model stage 1

v0 = [0;0];
b_glm = nan(NN,7);
yhat_glm = nan(NN,length(t));
yhat_slow = nan(NN,length(t));
llh_glm = nan(NN,1);

b_s1 = nan(NN,10);
yhat_s1 = nan(NN,length(t));
llh_s1 = nan(NN,1);
syn0_s1 = nan(NN,length(t));
syn_s1 = nan(NN,length(t));

for post= 1:NN
    
    
    y = CCG{post};
    
    % slow model
    b_glm(post,:) = glmfit(XX',y','poisson','constant','off');
    yhat_glm(post,:) = exp(b_glm(post,:)*XX);
    llh_glm(post) = nansum((y.*log( yhat_glm(post,:)+( yhat_glm(post,:)==0))- yhat_glm(post,:)));
    
    % full model Stage 1 (without latency constraints):
    
    b1 = random_parameter_s1(ACG,y,XX,t,b_glm(post,:),distance(post),eta,tau0,1); % 50 times random restart
    
    options=[];
    options.method = 'cg';
    options.MaxIter = 500;
    options.Display = 'off';
    [b,~,exitflag,output] = minFunc(@loss_excalpha,b1,options,ACG,XX',y',t',v0,distance(post),eta,tau0);
    [f,df,lam,syn] = loss_excalpha(b,ACG,XX',y',t',v0,distance(post),eta,tau0);
    b_s1(post,:) = b';
    yhat_s1(post,:) = lam';
    llh_s1(post) = nansum((y.*log(lam'+(lam'==0))-lam'));
    yhat_slow(post,:) = exp(b(1:size(XX,1))'*XX);
    syn_s1(post,:) = syn;
    
    p = exp(b((size(XX',2)+1):end));
    deltat = p(2);
    tau = p(3);
    t_syn = t;t_syn(t_syn<deltat)=deltat;
    syn0 = (t_syn-deltat)/tau.*exp(1-(t_syn-deltat)/tau);
    syn0 = syn0./( max(abs(syn0))+ (max(abs(syn0))==0));
    syn0_s1(post,:) = syn0;
    
end

llr = llh_s1 - llh_glm;

%% output

model_fits.glm.yhat = yhat_glm;
model_fits.glm.llh = llh_glm;
model_fits.glm.b = b_glm;

model_fits.inh = [];

model_fits.exc.yhat = yhat_s1;
model_fits.exc.yhat_slow = yhat_slow;
model_fits.exc.llh = llh_s1;
model_fits.exc.syn0 = syn0_s1; % without convolution
model_fits.exc.syn = syn_s1; % with convolution
model_fits.exc.b = b_s1;
model_fits.exc.eta = eta;
model_fits.info = 'stage1';


%% inhibitory model
% fprintf('inh model...')
% full model stage 1

% v0 = [0;0];
% b_glm = nan(NN,7);
% yhat_glm = nan(NN,length(t));
% yhat_slow = nan(NN,length(t));
% llh_glm = nan(NN,1);
%
% b_s1 = nan(NN,10);
% yhat_s1 = nan(NN,length(t));
% llh_s1 = nan(NN,1);
% syn0_s1 = nan(NN,length(t));
% syn_s1 = nan(NN,length(t));
%
% for post= 1:NN
%   
%   
%     y = CCG{post};
    
%    % full model Stage 1 (without latency constraints):
    
%     b1 = random_parameter_s1(ACG,y,XX,t,b_glm(post,:),distance(post),eta,tau0,1); % 50 times random restart
%     
%     options=[];
%     options.method = 'cg';
%     options.MaxIter = 500;
%     options.Display = 'off';
%     [b,~,exitflag,output] = minFunc(@loss_inhalpha,b1,options,ACG,XX',y',t',v0,distance(post),eta,tau0);
%     [f,df,lam,syn] = loss_inhalpha(b,ACG,XX',y',t',v0,distance(post),eta,tau0);
%     b_s1(post,:) = b';
%     yhat_s1(post,:) = lam';
%     llh_s1(post) = nansum((y.*log(lam'+(lam'==0))-lam'));
%     yhat_slow(post,:) = exp(b(1:size(XX,1))'*XX);
%     syn_s1(post,:) = syn;
%    
%     p = exp(b((size(XX',2)+1):end));
%     deltat = p(2);
%     tau = p(3);
%     t_syn = t;t_syn(t_syn<deltat)=deltat;
%     syn0 = (t_syn-deltat)/tau.*exp(1-(t_syn-deltat)/tau);
%     syn0 = syn0./( max(abs(syn0))+ (max(abs(syn0))==0));
%     syn0_s1(post,:) = syn0;
%     
% end
% 
% 
% llr = llh_s1 - llh_glm;
% 
% 
% %% output
% model_fits.inh.yhat = yhat_s1;
% model_fits.inh.yhat_slow = yhat_slow;
% model_fits.inh.llh = llh_s1;
% model_fits.inh.syn0 = syn0_s1; % without convolution
% model_fits.inh.syn = syn_s1; % with convolution
% model_fits.inh.b = b_s1;
% model_fits.inh.eta = eta;
% model_fits.info = 'stage1';





