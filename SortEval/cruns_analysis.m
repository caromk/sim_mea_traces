function [CT dist_from_probe_center] = cruns_analysis(filename,sampling_rate,len_seconds,varargin)

% get options 
opt = Opt(varargin);

% defaults
default_spike_sort = 0;

% can run with a shorter length
spike_sort = opt.get('SpikeSort',default_spike_sort);

CT.sampling_rate = sampling_rate;
CT.loc_neuron = csvread(sprintf('%s%s',filename,'_neuron_loc.txt'));

% spike times
spikes = csvread(sprintf('%s%s',filename,'_neuron.txt'));
CT.n_neurons = size(spikes,1);
for i_neuron = 1:CT.n_neurons
    CT.spike_times{i_neuron} = spikes(i_neuron,spikes(i_neuron,:)>0)';
end


CT.traces = csvread(sprintf('%s%s',filename,'_detector.txt'));
CT.loc_electrode = csvread(sprintf('%s%s',filename,'_detector_loc.txt'));
CT.n_electrode = size(CT.loc_electrode,1);

if size(CT.traces,2) == 1
    CT.traces = reshape(CT.traces,CT.sampling_rate*len_seconds,[])';
end

% dist_from_soma
CT.dist_from_soma = zeros(CT.n_electrode,CT.n_neurons);
for i_neuron = 1:CT.n_neurons
    for i_electrode = 1:CT.n_electrode
        CT.dist_from_soma(i_electrode,i_neuron) = norm(CT.loc_neuron(i_neuron,:) - CT.loc_electrode(i_electrode,:));
    end
end

dist_from_probe_center = sum(CT.loc_neuron.^2,2).^(1/2);

if spike_sort
    [E comps all_thres_spikes spike_comps wts] = SortEval_cruns(CT);
    save(filename,'CT','dist_from_probe_center','E','comps','all_thres_spikes','spike_comps','wts','-v7.3');
else
    save(filename,'CT','dist_from_probe_center','-v7.3');
end
