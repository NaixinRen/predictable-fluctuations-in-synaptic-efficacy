function [efficacy,effstd] = efficacy_estimate(chunkdata,GLM_mdlfit_chunk)

numspk_pre = chunkdata.numspk_pre;
nchunk = length(numspk_pre);

yhat = GLM_mdlfit_chunk.exc.yhat;
yhat_slow = GLM_mdlfit_chunk.exc.yhat_slow;
efficacy = sum(yhat - yhat_slow,2)./numspk_pre';


effstd = nan(nchunk,1);
for i = 1:nchunk
    if isnan(sum(GLM_mdlfit_chunk.exc.b(i,:)))
        continue
    end
    bsamp = mvnrnd(GLM_mdlfit_chunk.exc.b(i,:),GLM_mdlfit_chunk.exc.covb{i},100);
    yhatsamp = exp(bsamp*GLM_mdlfit_chunk.exc.basis');
    yslowsamp = exp(bsamp(:,1)*GLM_mdlfit_chunk.exc.basis(:,1)');
    effsamp = sum(yhatsamp-yslowsamp,2)./numspk_pre(i);
    effstd(i) = std(effsamp);
end