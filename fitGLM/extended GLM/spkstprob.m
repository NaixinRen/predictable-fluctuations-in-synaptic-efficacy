
function [ignore_spkst,SortingProbIndex] = spkstprob(CCG,interval,binsize)
n = size(CCG,2);
SortingProbIndex = nan(n,n);
%SortingProbIndex_c = nan(n,n);
ignore_spkst = zeros(n,n);
c = interval/(2*binsize)+1;
nbin = min(floor(1.5/binsize),5);
for pre = 1:n
    for post = 1:n
        if isempty(CCG{pre,post}) 
            continue
        end       
        mu(1) = mean(CCG{pre,post}(c-(nbin-1)/2:c+(nbin-1)/2)); %center 1.5 ms
        mu(2) = mean(CCG{pre,post}(c-nbin-(nbin-1)/2:c-(nbin-1)/2-1)); %left 1.5 ms
        mu(3) = mean(CCG{pre,post}(c+(nbin-1)/2+1:c+nbin+(nbin-1)/2)); %right 1.5 ms
        SortingProbIndex(pre,post) = nanmin( (mu(1)-mu(2))/mu(2), (mu(1)-mu(3))/mu(3) );
        %SortingProbIndex(post,pre) = SortingProbIndex(pre,post);
        %SortingProbIndex_c(pre,post) = CCG{pre,post}(c)/min(CCG{pre,post}(c-1),CCG{pre,post}(c+1));
    end
end
ignore_spkst(SortingProbIndex >= 1) = 1;
ignore_spkst(isnan(SortingProbIndex)) = 1;
