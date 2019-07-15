cd  /home/caromk/simdata/rand_med1/
filename = 'rand_med1M';
load(filename);

addpath('~/neuranalysis/SpikeSort/')
addpath('~/neuranalysis/Opt/')
addpath('~/neuranalysis/SimTraces/SortEval/')

for count = [50 1700:100:2000 123 100:100:400]
    CT = subsetT(CT,1:count);
    save('-v7.3',sprintf('%s_%d.mat',filename,count),'CT','dist_from_probe_center');
    [E comps all_thres_spikes spike_comps wts] = SortEval_cruns(CT);
    save('-v7.3',sprintf('%s_%d.mat',filename,count),'E','comps','all_thres_spikes','spike_comps','wts','-append');
    clear all
filename = 'rand_med1M';
    load(filename);
end
