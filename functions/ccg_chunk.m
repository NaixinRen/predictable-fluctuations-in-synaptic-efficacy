function [ccg,fr_pre,fr_post,numspk_pre,numspk_post,accg] = ccg_chunk(spk_pre,spk_post,info)

chunktime = info.chunktime;
overlap = info.chunktime_overlap;
%nchunk = floor(min(spk_pre(end),spk_post(end))/(60*chunktime));
nchunk = floor( (min(max(spk_pre),max(spk_post))-60*chunktime)/(60*(chunktime-overlap)) )+1;
spk_prechunk = cell(1,nchunk);spk_postchunk = cell(1,nchunk);
fr_pre = zeros(1,nchunk);fr_post = zeros(1,nchunk);
numspk_pre = zeros(1,nchunk);numspk_post = zeros(1,nchunk);
for n = 1:nchunk
%     t1 = .001+(n-1)*(60*chunktime);
%     t2 = n*(60*chunktime);

    t1 = .001+(n-1)*(chunktime-overlap)*60;
    t2 = (chunktime + (n-1)*(chunktime-overlap))*60;
    
    temp = spk_pre>t1 & spk_pre<t2;
    fr_pre(n) = sum(temp)/(60*chunktime);
    numspk_pre(n) = sum(temp);
    spk_prechunk{n} = spk_pre(temp);
    
    temp = spk_post>t1 & spk_post<t2;
    fr_post(n) = sum(temp)/(60*chunktime);
    numspk_post(n) = sum(temp);
    spk_postchunk{n} = spk_post(temp);
end

interval = info.interval;
binsize = info.binsize;

ccg = cell(1,nchunk);
for n = 1:nchunk
    if isempty(spk_prechunk{n}) || isempty(spk_postchunk{n})
        continue
    end
    [d,~] = corr_fast_v3(spk_prechunk{n}*1000,spk_postchunk{n}*1000,-interval/2,interval/2,interval/binsize+2);
    [ad,~] = corr_fast_v3(spk_prechunk{n}*1000,spk_prechunk{n}*1000,-interval,interval,2*interval/binsize+2);
    
    if isempty(d)
        continue
    end
    ccg{n}(1,:) = d(1:interval/binsize+1)';
    accg{n}(1,:) = ad(1:end-1)';
    
    
    
end