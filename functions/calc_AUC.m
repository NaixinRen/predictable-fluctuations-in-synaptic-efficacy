function [AUC,AUC_spon] = calc_AUC(s_pre, s_post, mdlfit, fit)
prespklist = find(s_pre);
yhat = mdlfit.yhat;
label_true = nan(length(prespklist),1);
score = nan(length(prespklist),1);
k = find(fit.Fc>.5) - ceil(length(fit.Fc)/2);
for i = 1:length(prespklist)
    range = prespklist(i)+k;
    range = range(range<length(s_post));
    y_i = s_post(range);
    yhat_i = yhat(range);
    label_true(i) = sum(y_i)>0;
    score(i) = 1-prod(poisscdf(0,yhat_i)); 
end

[X,Y,T,AUC] = perfcurve(label_true,score,1);
dt=0.001;
prelist_spon = find(prespklist<=6190/dt & prespklist>4400/dt); % spontaneous activity

if sum(label_true(prelist_spon))==0 || sum(label_true(prelist_spon))== length(label_true(prelist_spon))
    AUC_spon = nan;
else
    [X,Y,T,AUC_spon] = perfcurve(label_true(prelist_spon),score(prelist_spon),1);
end