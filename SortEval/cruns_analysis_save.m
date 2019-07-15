function cruns_analysis_save(filename,sampling_rate,length)

if isa(sampling_rate,'char')
    sampling_rate = str2double(sampling_rate);
end

if isa(length,'char')
    length = str2double(length);
end 

disp(sprintf('%s Loading files.',datestr(now)));

CT.sampling_rate = sampling_rate;

CT.traces = csvread(sprintf('%s%s',filename,'_detector.txt'));
CT.loc_neuron = csvread(sprintf('%s%s',filename,'_neuron_loc.txt'));
CT.loc_electrode = csvread(sprintf('%s%s',filename,'_detector_loc.txt'));
CT.n_electrode = size(CT.loc_electrode,1);

if size(CT.traces,2) == 1
    CT.traces = reshape(CT.traces,sampling_rate*length,[])';
end

% spike times
spikes = csvread(sprintf('%s%s',filename,'_neuron.txt'));
CT.n_neurons = size(spikes,1);
for i_neuron = 1:CT.n_neurons
    CT.spike_times{i_neuron} = spikes(i_neuron,spikes(i_neuron,:)>0)';
end

% dist_from_soma
CT.dist_from_soma = zeros(CT.n_electrode,CT.n_neurons);
for i_neuron = 1:CT.n_neurons
    for i_electrode = 1:CT.n_electrode
        CT.dist_from_soma(i_electrode,i_neuron) = norm(CT.loc_neuron(i_neuron,:) - CT.loc_electrode(i_electrode,:));
    end
end

dist_from_probe_center = sum(CT.loc_neuron.^2,2).^(1/2);
save('-v7.3',sprintf('%s.mat',filename),'CT','dist_from_probe_center');

disp(sprintf('%s Mat file saved.',datestr(now)));

[E comps all_thres_spikes spike_comps wts] = SortEval_cruns(CT);

save('-v7.3',sprintf('%s.mat',filename),'E','comps','all_thres_spikes','spike_comps','wts','-append');
