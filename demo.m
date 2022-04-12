clear;clc;close all
addpath(genpath('fitGLM'))
addpath(genpath('functions'))
addpath(genpath('GBLM'))

load('spikes_from_715093703.mat')
%% generate CCG & ACG& fit GLM (estimate efficacy)
dt=0.001;
Info.interval = 50;Info.binsize = 1; %ms
% CCG
[tmp,~] = corr_fast_v3(Tpre/dt,Tpost/dt,-Info.interval/2,Info.interval/2,(Info.interval/Info.binsize)+2);
CCG{1} = tmp(1:end-1)';

% double sized ACG
[tmp,~] = corr_fast_v3(Tpre/dt,Tpre/dt,-Info.interval,Info.interval,2*(Info.interval/Info.binsize)+2);
ACG{1} = tmp(1:end-1)';


figure(1),
subplot(2,2,1)
t = linspace(-25,25,length(CCG{1}));
bar(t,CCG{1},'k')
hold on
plot([sum(xlim)/2 sum(xlim)/2],ylim,'k--')
box off
title('CCG')
xlabel('interval [ms]')
subplot(2,2,2)
bar(linspace(-50,50,length(ACG{1})),ACG{1},'k')
hold on
plot([sum(xlim)/2 sum(xlim)/2],ylim,'k--')
box off
title('ACG')
xlabel('interval [ms]')

% divide spike trains to 5-min windows
Info.chunktime = 5; % mins
Info.chunktime_overlap = Info.chunktime-1;
[CCG_chunk,fr_pre,fr_post,numspk_pre,numspk_post,ACG_chunk] = ccg_chunk(Tpre,Tpost,Info);
chunkdata.CCG = CCG_chunk;
chunkdata.fr_pre = fr_pre;
chunkdata.fr_post = fr_post;
chunkdata.numspk_pre = numspk_pre;
chunkdata.numspk_post = numspk_post;
chunkdata.info = Info;
% fit GLM to estimated efficacy
load('learned_basis.mat')
hyperparameter.tau0 = 0.8; %ms
hyperparameter.eta_w = 0;
hyperparameter.eta_tau = 10;
GLM_mdlfit = extendedGLM_single(CCG,ACG{1}, X, hyperparameter,Info);
GLM_mdlfit_chunk = extendedGLM_chunk(CCG_chunk, X, hyperparameter,Info,GLM_mdlfit);
[GLM_efficacy,GLM_effstd] = efficacy_estimate(chunkdata,GLM_mdlfit_chunk);


subplot(2,2,1)
hold on
plot(t,GLM_mdlfit.exc.yhat,'linewidth',1)
hold on
plot(t,GLM_mdlfit.exc.yhat_slow,'linewidth',1)

subplot(2,2,3)
Fc = GLM_mdlfit.exc.syn0;
syn_delay = exp(GLM_mdlfit.exc.b(end-1));
plot(t,Fc,'linewidth',1)
title('Coupling filter')
box off

subplot(2,2,4)
nknots=4;
Xslow = getCubicBSplineBasis(linspace(-2/(nknots-1),1+2/(nknots-1),151),4,false);
Xslow = Xslow(:,2:end); % get rid of constant term
plot(-floor(size(Xslow,1)/2):floor(size(Xslow,1)/2),Xslow)
title('Slow effect filter')
box off
set(gcf,'position',[200,200,500,400])

figure(2),clf;
subplot(2,1,1)
nchunk = length(fr_pre);chunkwidth = Info.chunktime-Info.chunktime_overlap;
t = (0:nchunk-1)*chunkwidth+Info.chunktime/2;
plot(t,fr_pre,'k','Linewidth',2)
box off
ylabel(['Firing rate [Hz]'])
subplot(2,1,2)
fill([t,flip(t)],[GLM_efficacy'-GLM_effstd',flip(GLM_efficacy'+GLM_effstd')],[0, 0.4470, 0.7410],'FaceAlpha',.3,'EdgeColor','none')
hold on
p1 = plot(t,GLM_efficacy,'Color',[0, 0.4470, 0.7410],'linewidth',2);
box off
xlabel('Time [min]')
ylabel('Efficacy')

set(gcf,'position',[700,200,500,400])
%% ISI - efficacy
dt=0.001;
ISI = diff(Tpre);
edges= linspace(min(log(ISI)),max(log(ISI)),21);
isidata.ccg_isi = [];isi_m=[];isidata.numspk_pre = [];
for i = 1:20
    list = find(log(ISI)>=edges(i) & log(ISI) <edges(i+1));
    Tpre_isi = Tpre(list+1);
    [d,~] = corr_fast_v3(Tpre_isi/dt,Tpost/dt,-Info.interval/2,Info.interval/2,Info.interval/Info.binsize+2);
    if isempty(d)
        continue
    end
    d = reshape(d,[Info.interval/Info.binsize+2,1]);
    isidata.ccg_isi{i} = d(1:Info.interval/Info.binsize+1)';
    isidata.numspk_pre(1,i) = length(list);
    isi_m(i) = mean([exp(edges(i)),exp(edges(i+1))]);
end

GLM_mdlfit = extendedGLM_single(CCG,ACG{1}, X, hyperparameter,Info);
GLM_mdlfit_chunk = extendedGLM_chunk(isidata.ccg_isi, X, hyperparameter,Info,GLM_mdlfit);
[efficacy_isi,effstd] = efficacy_estimate(isidata,GLM_mdlfit_chunk);
efficacy_isi( isidata.numspk_pre<50)=nan;
efficacy_isi(efficacy_isi<0)=0;

figure(3),clf
subplot(2,1,1)
[~,edges] = histcounts(log10(ISI/dt),500);
 histogram(ISI/dt,10.^edges,'FaceColor','k','FaceAlpha',.5,'EdgeColor','none')
 set(gca, 'xscale','log')
xticks([1 10 100 1000 10000])
xticklabels({'10^0','10^1','10^2','10^3','10^4'})
set(gca, 'XScale', 'log')
xlim([1 10000])
set(gca,'linewidth',1)
set(gca,'FontSize',12)
ax1 = gca;
ax1.YAxis.Visible = 'off';
box off
subplot(2,1,2)
scatter(isi_m/dt,efficacy_isi,'MarkerFaceColor','k','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
xticks([1 10 100 1000 10000])
xticklabels({'10^0','10^1','10^2','10^3','10^4'})
set(gca, 'XScale', 'log')
xlim([1 10000])
box off
set(gca,'linewidth',1)
set(gca,'FontSize',12)
xlabel('Presynaptic ISI [ms]')
ylabel('Synaptic efficacy')
grid on
set(gcf,'position',[100,200,300,300])
%% run GBLM
population{1}=Tpre;
population{2}=Tpost;
S = double(getSpkMat(population,dt,[],1));
range = 1:length(S);
s_pre = S(1,range);s_post = S(2,range);

fit.tau=0.2;
fit.Nq = 8;
fit.hist_beta = 0;
fit.Fc = Fc;
fit.Fsmooth = Xslow';
fit.F = eye(1);
fit.Q = eye(1)*1e-3;
toleranceValue = 5;
% GBLM
[theta,mdlfit] = ...
    stp_gblm_adaptive_v5(s_pre,s_post,fit,'toleranceValue',toleranceValue);
% model with static synaptic effect
[theta_s,mdlfit_s] = ...
    static_glm_adaptive(s_pre,s_post,fit,'toleranceValue',toleranceValue);

[AUC,AUC_spon] = calc_AUC(s_pre, s_post, mdlfit, fit);

%%%% GBLM efficacy %%%%
nchunk = length(CCG_chunk);
chunkwidth = Info.chunktime-Info.chunktime_overlap;
yhat_pre = mdlfit.yhat-mdlfit.yhat0;
yhat_pre_static = mdlfit_s.yhat-mdlfit_s.yhat0;

GBLM_eff = nan(nchunk,1);GBLM_eff_static = nan(nchunk,1);
for i = 1:nchunk
    t1 = 1+(i-1)*chunkwidth*60/dt;
    t2 = (Info.chunktime + (i-1)*chunkwidth)*60/dt;
    if min(t1,t2)>length(yhat_pre)
        break;
    end
    GBLM_eff(i) = sum(yhat_pre(t1:t2))/numspk_pre(i);
    GBLM_eff_static(i) = sum(yhat_pre_static(t1:t2))/numspk_pre(i);
end
figure(2),
subplot(2,1,2)
p2 = plot(t,GBLM_efficacy,'color',[0.8500, 0.3250, 0.0980],'linewidth',2);
l1 = legend([p1,p2],{'observed','GBLM estimated'});
set(l1,'Box','off')
%% modi function
modif = theta.modif_fxn*theta.w;
se_modif = theta.se_modif_fxn*theta.w;
isi = 1:length(modif);
figure(4),clf
subplot(2,1,1)
dt=.001;
ISI = diff(Tlist{siglist.pre(syn)});
edges = logspace(0,4,501);
histogram(ISI/dt,edges,'FaceColor','k','FaceAlpha',.5,'EdgeColor','none')
xticks([1 10 100 1000 1e4])
xticklabels({'10^0','10^1','10^2','10^3','10^4'})
set(gca, 'XScale', 'log')
xlim([1 10000])
box off
set(gca,'linewidth',1)
set(gca,'FontSize',12)
ax1 = gca;
ax1.YAxis.Visible = 'off';
subplot(2,1,2)
fill([isi,flip(isi)],[modif-se_modif,flip(modif+se_modif)],[0.8500, 0.3250, 0.0980]	,'FaceAlpha',.3,'EdgeColor','none')
hold on
plot(isi,modif,'Color',[0.8500, 0.3250, 0.0980]	,'linewidth',2)
hold on
xlim([1 10000])
plot(xlim,[0,0],'k--')
xticks([1 10 100 1000 1e4])
xticklabels({'10^0','10^1','10^2','10^3','10^4'})
grid on
% ylim([-.5,.5])
xlabel('ISI[ms]')
box off
set(gca, 'XScale', 'log')
set(gca,'linewidth',1)
set(gca,'FontSize',12)
set(gcf,'position',[400,200,300,300])
%% shuffling

params = struct();
params.trange = [-.02 .02];
params.tn = 202;
params.dt = range(params.trange)/(params.tn-2);
[c,dT] = corr_fast_v3(Tpre,Tpost,params.trange(1)-params.dt/2,params.trange(2)+params.dt/2,params.tn);

R_pre = nan(100,1);R_post = nan(100,1);CV_shuff = nan(100,1);
for shuff_n = 1:10
    shuff_n
    keepidx=zeros(size(dT,1),1);
    ua = unique(dT(:,3));
    for i=1:length(ua)
        r = find(dT(:,3)==ua(i));
        keepidx(datasample(r,1))=1;
    end
    n = length(ua);
    ShuffPost = Tpost;
    ShuffPost = Tpost(setdiff(1:length(ShuffPost),unique(dT(:,3)))); % keep spikes that aren't in the xcorr
    presample = datasample(Tpre,n);
    intsample = datasample(dT(keepidx>0,1),n);
    ShuffPost = sort([ShuffPost, presample+intsample']);
    Tpost_s=ShuffPost;

    [CCG_chunk,fr_pre,fr_post,numspk_pre,numspk_post] = ccg_chunk(Tpre,Tpost_s,Info);
    chunkdata = struct();
    chunkdata.CCG = CCG_chunk;chunkdata.info = Info;
    chunkdata.fr_pre = fr_pre;chunkdata.fr_post = fr_post;
    chunkdata.numspk_pre = numspk_pre;chunkdata.numspk_post = numspk_post;

    GLM_mdlfit_chunk = extendedGLM_chunk(CCG_chunk, X, hyperparameter,Info,GLM_mdlfit);
    [efficacy,effstd] = efficacy_estimate(chunkdata,GLM_mdlfit_chunk);
    efficacy(efficacy<0)=0;
    [R_pre(shuff_n,:),~] = corr(efficacy(~isnan(efficacy)),fr_pre(~isnan(efficacy))', 'Type','Spearman');
    [R_post(shuff_n,:),~] = corr(efficacy(~isnan(efficacy)),fr_post(~isnan(efficacy))', 'Type','Spearman');
    CV_shuff(shuff_n,:) = nanstd(efficacy)/nanmean(efficacy);
end
figure(5),clf
histogram(R_pre,100,'FaceColor','k','EdgeColor','none','FaceAlpha',.5)
xlabel('Spearman correlation')
hold on
[R_true_pre,~] = corr(GLM_efficacy(~isnan(GLM_efficacy)),fr_pre(~isnan(GLM_efficacy))', 'Type','Spearman');
plot([R_true_pre,R_true_pre],ylim,'linewidth',2)
xlim([-1,1])
box off
set(gca,'linewidth',1)
set(gca,'FontSize',12)
set(gcf,'position',[1200,400,200,200])