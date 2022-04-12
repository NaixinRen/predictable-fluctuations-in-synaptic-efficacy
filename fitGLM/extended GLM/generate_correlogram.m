function [CCG,ACCG, ignore] = generate_correlogram(Tlist,sr,interval,binsize)


% This function generates cross-correlograms from spike trains

% input arguments:
% -spikes: spike times
% -sr: sampling rate of the spike data (per ms) 
% -location: 2D location of the neurons(soma)
% -interval: interval of the correlogram. eg. interval = 50 means the interval of the correlogram is [-25,25] ms
% -ignore_index: 1 or 0, when ignore_index == 1,ignore low firing rate neurons (# of spikes <1000), sparse correlograms 
% (# of spikes in a correlogram <50), and neuron pairs with possbile spike sorting problems. Defalt value is 1.


NN = length(Tlist); 

%Tlist=cell(0);
for i=1:NN
    %add sub-resolution noise
    Tlist{i}=Tlist{i}+(rand(size(Tlist{i}))-.5)*(1/sr);
end


CCG = cell(NN);
for i=1:NN
    fprintf('neuron %i \n',i)
    if length(Tlist{i}) <= 1000
        continue
    end
%     [tmp,~] = corr_fast_v3(Tlist{i},Tlist{i},-interval/2,interval/2,2*(interval/binsize)+2); - wrong code
    [tmp,~] = corr_fast_v3(Tlist{i},Tlist{i},-interval,interval,2*(interval/binsize)+2);
    ACCG{i} = tmp(1:end-1)'; % last bin of histc contains data == the last bin
    for j=(i+1):NN
        if length(Tlist{j}) <= 1000
            continue
        end
        [tmp,~] = corr_fast_v3(Tlist{i},Tlist{j},-interval/2,interval/2,interval/binsize+2);
        CCG{i,j} = tmp(1:interval/binsize+1)'; % last bin of histc contains data == the last bin
        CCG{j,i} = flipud(tmp(1:interval/binsize+1))';
    end
end

ignore = eye(NN,NN);
ignore(cellfun(@sum,CCG) < 50) = 1;
