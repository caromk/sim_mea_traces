function [E comps varargout] = SortEval(T,varargin)

default_alt_length = size(T.traces,2);
default_comps = [];
default_verbose = 0;
default_peak_variation_ratio = 0.2;
default_zero_interval_separation_min = 0.2;
default_compare_interval_sec = 0.8e-3;
default_all_thres_spikes = {};

% get options 
opt = Opt(varargin);

% can run with a shorter length
alt_length = opt.get('AltLength',default_alt_length);
% init optional inputs / defaults
% verbose prints out details and displays figures
verbose = opt.get('Verbose',default_verbose);
% ratio for allowed variation in peak height, in ratio of peak height
peak_variation_interval = opt.get('PeakVariationRatio',default_peak_variation_ratio);
% minimum separation distance (in max_peak scale) between nonzero intervals
% for the last nonzero interval to be the robust peaks
zero_interval_separation_min = opt.get('ZeroIntervalSeparationMin',default_zero_interval_separation_min);
% compare interval is the window used to compare spikes and spike trains to each other to see if they're the same 
compare_interval_sec = opt.get('CompareIntervalSec',default_compare_interval_sec);
% thres spikes contains date to decrease how much needs to be rerun
all_thres_spikes = opt.get('AllThresSpikes',default_all_thres_spikes);
% comps = components returned from ICA, added as an option for cases in
% which the algorithm needs to be run multiple times, makes it faster and
% more consistent
comps = opt.get('Comps',default_comps);
if isempty(comps)
    fprintf('Running comps.\n')
end

tic
[E.spike_trains,spike_comps,comps,E.n_dups,E.best_spike_index,E.ind_spiking_comps,all_thres_spikes] = RobustSpikeSort(T.traces(:,1:alt_length),T.sampling_rate,'Verbose',verbose,'PeakVariationRatio',peak_variation_interval,'ZeroIntervalSeparationMin',zero_interval_separation_min,'Comps',comps,'AllThresSpikes',all_thres_spikes);
E.run_time = toc;
E.peak_variation_interval = peak_variation_interval;
E.zero_interval_separation_min = zero_interval_separation_min;

E.n_spiking_comp = size(E.ind_spiking_comps,2);
E.matches = [];
E.best_match = zeros(E.n_spiking_comp,6);
E.n_spiking_comp_accurate = 0;
E.rate_false_pos_spikes = zeros(E.n_spiking_comp,1);
E.rate_false_neg_spikes = zeros(E.n_spiking_comp,1);
E.n_false_pos_spikes = zeros(E.n_spiking_comp,1);
E.n_false_neg_spikes = zeros(E.n_spiking_comp,1);

disp(sprintf('%s %d spiking components',datestr(now),E.n_spiking_comp))
for i_spiking_comp = 1:E.n_spiking_comp
    %% format of curr_matches/best_match: i_neuron n_putative_spikes match_rate_puta match_rate_real match_mean_diff match_mean_std
    curr_matches = FindBestSpikeTimesMatch(E.spike_trains(i_spiking_comp,:),T.spike_times,compare_interval_sec*T.sampling_rate);
    E.matches{i_spiking_comp} = curr_matches;
    if ~isempty(curr_matches)
        if size(find(curr_matches(curr_matches(:,3) == 1,4) == 1),1) == 1
            E.n_spiking_comp_accurate = E.n_spiking_comp_accurate + 1;
            E.best_match(i_spiking_comp,:) = curr_matches(curr_matches(:,3) == 1 & curr_matches(:,4) == 1,:);
        elseif size(find(curr_matches(:,3) == 1),1) == 1
            E.best_match(i_spiking_comp,:) = curr_matches(curr_matches(:,3) == 1,:);
            n_spikes_orig = E.best_match(i_spiking_comp,2)/E.best_match(i_spiking_comp,4);
            n_spikes_comp = E.best_match(i_spiking_comp,2);
            spike_index = T.spike_times{E.best_match(i_spiking_comp,1)};
            if n_spikes_orig - n_spikes_comp == 2
                if spike_index(1) - compare_interval_sec*T.sampling_rate - E.best_match(i_spiking_comp,5) - 2*E.best_match(i_spiking_comp,6) < 0 && spike_index(end) + compare_interval_sec*T.sampling_rate + E.best_match(i_spiking_comp,5) + 2*E.best_match(i_spiking_comp,6) > size(E.spike_trains(i_spiking_comp,:),2)
                    n_spiking_comp_accurate = n_spiking_comp_accurate + 1;
                end
            end
            if n_spikes_orig - n_spikes_comp == 1
                if spike_index(1) - compare_interval_sec*T.sampling_rate - E.best_match(i_spiking_comp,5) - 2*E.best_match(i_spiking_comp,6) < 0 || spike_index(end) + compare_interval_sec*T.sampling_rate + E.best_match(i_spiking_comp,5) + 2*E.best_match(i_spiking_comp,6) > size(E.spike_trains(i_spiking_comp,:),2)
                    E.n_spiking_comp_accurate = E.n_spiking_comp_accurate + 1;
                end
            end
        elseif size(find(curr_matches(:,3) == 1),1) > 1
            curr_matches_puta_exact = curr_matches(curr_matches(:,3) == 1,:);
            [val_max_real i_max_real] = max(curr_matches_puta_exact(:,4));
            E.best_match(i_spiking_comp,:) = curr_matches_puta_exact(i_max_real,:);
        else
            [val_max_puta i_max_puta] = max(curr_matches(:,3));
            E.best_match(i_spiking_comp,:) = curr_matches(i_max_puta,:);
        end
    end
    E.rate_false_pos_spikes(i_spiking_comp) = 1 - E.best_match(i_spiking_comp,3);
    E.rate_false_neg_spikes(i_spiking_comp) = 1 - E.best_match(i_spiking_comp,4);
    E.n_false_pos_spikes(i_spiking_comp) = (1 - E.best_match(i_spiking_comp,3))*E.best_match(i_spiking_comp,2);
    E.n_false_neg_spikes(i_spiking_comp) = E.best_match(i_spiking_comp,2)/E.best_match(i_spiking_comp,4) - E.best_match(i_spiking_comp,2);
end
disp(sprintf('%s %d matches found!',datestr(now),E.n_spiking_comp_accurate))

E.n_neurons_50microV = sum(sum(T.calc_snr_multiplier(T.dist_from_soma)>50)>=1);
E.min_dist_from_pad = zeros(E.n_spiking_comp,1);
E.min_dist_from_pad(E.best_match(:,1)>=1) = min(T.dist_from_soma(:,E.best_match(E.best_match(:,1)>=1,1)))';

y_max = max(T.loc_electrode(:,2));
E.min_dist_from_probe_center_line = zeros(E.n_spiking_comp,1);
goodmatches = sum(E.best_match(:,1)>=1);
E.min_dist_from_probe_center_line(E.best_match(:,1)>=1) = sqrt(T.loc_neuron(E.best_match(E.best_match(:,1)>=1,1),1).^2 + max(zeros(goodmatches,1),abs(T.loc_neuron(E.best_match(E.best_match(:,1)>=1,1),2))-y_max).^2 + T.loc_neuron(E.best_match(E.best_match(:,1)>=1,1),3).^2);

varargout{1} = all_thres_spikes;