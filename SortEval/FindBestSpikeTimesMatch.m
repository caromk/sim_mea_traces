% input:
%   test_spike_times - single list of spike times, either as 0s and 1s for 
%     the length of the trace with a 1 at the spike peak or indices eg find
%   match_spike_times - set of lists of spike times, either as a cell array
%     if given as indices or a matrix if given as 0s and 1s
%   compare_interval - maximum number of timepoints the peaks can be apart
%     to be considered the same spike
% optional input: (specify with flag followed by value)
%   'AltLength' - match on a shorter number of time points
%   'MaxSpikeDiffCompare' - to cut processing time, can limit the difference
%     in number of spikes between any evaluated pair of spike times

% todo: change from n_neurons to n_matches or something for clarity for
% users who are using for other uses (no clarity in this sentence tho)

function matches = FindBestSpikeTimesMatch(test_spike_times,match_spike_times,compare_interval,varargin)

default_alt_length = 0;
default_max_spike_diff_compare = 5;
default_min_match_rate_to_report = 0.1;

% get options 
opt = Opt(varargin);
% set options
alt_length = opt.get('AltLength',default_alt_length);
max_spike_diff_compare = opt.get('MaxSpikeDiffCompare',default_max_spike_diff_compare);
min_match_rate_to_report = opt.get('MinMatchRateToReport',default_min_match_rate_to_report);

% check if test_spike_times is indices or 0s and 1s
uniq_vals = unique(test_spike_times);
if (size(uniq_vals,2)==2 && ismember(0,uniq_vals) && ismember(1,uniq_vals)) || (size(uniq_vals,2)==1 && size(uniq_vals,1)==1 && (uniq_vals==1 || uniq_vals==0))
    test_spike_times_index = find(test_spike_times);
else
    test_spike_times_index = test_spike_times;
end

% make sure 
if size(test_spike_times_index,1)~=1
    test_spike_times_index = test_spike_times_index';
end

% number of test spikes
n_test_spikes = size(test_spike_times_index,2);

% check if match_spike_times is in cell format, if not make it so
if ~iscell(match_spike_times)
    n_neurons = size(match_spike_times,1);
    old_match_spike_times = match_spike_times;
    match_spike_times = {};
    for i_neuron = 1:n_neurons
        match_spike_times{i_neuron} = find(old_match_spike_times(i_neuron,:))';
    end
end

% number of neurons
n_neurons = size(match_spike_times,2);

if alt_length > 0
    for i_neuron = 1:n_neurons
        match_spike_times{i_neuron} = match_spike_times{i_neuron}(match_spike_times{i_neuron}<=alt_length);
    end
end

matches = zeros(n_neurons,6);
i_match = 0;

if ~isempty(test_spike_times)
    for i_neuron = 1:n_neurons
        if abs(n_test_spikes-size(match_spike_times{i_neuron},1)) < max_spike_diff_compare
            [ab_pair ba_pair a_unpair b_unpair] = CompareTwoSpikeTimes(test_spike_times_index,match_spike_times{i_neuron},compare_interval);
            n_match_spikes = length(match_spike_times{i_neuron});
            pairs = [ab_pair ba_pair];
            match_rate_test = size(pairs,2)/n_test_spikes;
            match_rate_match = size(pairs,2)/n_match_spikes;
            if match_rate_test >= min_match_rate_to_report
                i_match = i_match + 1;
                match_mean_diff = mean(pairs(1,:)-pairs(2,:));
                match_mean_std = std(pairs(1,:)-pairs(2,:));
                matches(i_match,:) = [i_neuron n_test_spikes match_rate_test match_rate_match match_mean_diff match_mean_std];
            end
        end
    end
    
end

if i_match == 0
    matches = [];
else
    matches = matches(1:i_match,:);
end