function [matches,best_match,n_spiking_comp_accurate,rate_false_pos_spikes,rate_false_neg_spikes,n_false_pos_spikes,n_false_neg_spikes,accurate] = EvalSorted(T,spike_trains,compare_interval_sec,varargin)

default_alt_length = size(spike_trains,2);

% get options 
opt = Opt(varargin);

% can run with a shorter length
alt_length = opt.get('AltLength',default_alt_length);

n_spiking_comp = size(spike_trains,1);
matches = [];
best_match = zeros(n_spiking_comp,6);
n_spiking_comp_accurate = 0;
rate_false_pos_spikes = zeros(n_spiking_comp,1);
rate_false_neg_spikes = zeros(n_spiking_comp,1);
n_false_pos_spikes = zeros(n_spiking_comp,1);
n_false_neg_spikes = zeros(n_spiking_comp,1);
accurate = zeros(n_spiking_comp,1);

for i_spiking_comp = 1:n_spiking_comp
    % format of curr_matches/best_match: 1-i_neuron 2-n_putative_spikes 3-match_rate_puta 4-match_rate_real 5-match_mean_diff 6-match_mean_std
    curr_matches = FindBestSpikeTimesMatch(spike_trains(i_spiking_comp,:),T.spike_times,compare_interval_sec*T.sampling_rate,'AltLength',alt_length);
    matches{i_spiking_comp} = curr_matches;
    if ~isempty(curr_matches)
        % one perfect match
        if size(find(curr_matches(curr_matches(:,3) == 1,4) == 1),1) == 1
            n_spiking_comp_accurate = n_spiking_comp_accurate + 1;
	    accurate(i_spiking_comp) = 1;
            best_match(i_spiking_comp,:) = curr_matches(curr_matches(:,3) == 1 & curr_matches(:,4) == 1,:);
        % every sorted spike is found in the real data, but real data
        % spikes are unaccounted for
        elseif size(find(curr_matches(:,3) == 1),1) == 1
            best_match(i_spiking_comp,:) = curr_matches(curr_matches(:,3) == 1,:);
            n_spikes_orig = round(best_match(i_spiking_comp,2)/best_match(i_spiking_comp,4));
            n_spikes_comp = best_match(i_spiking_comp,2);
            spike_index = T.spike_times{best_match(i_spiking_comp,1)};
            if n_spikes_orig - n_spikes_comp == 2
                if spike_index(1) - compare_interval_sec*T.sampling_rate - best_match(i_spiking_comp,5) - 2*best_match(i_spiking_comp,6) < 0 && spike_index(end) + compare_interval_sec*T.sampling_rate + best_match(i_spiking_comp,5) + 2*best_match(i_spiking_comp,6) > size(spike_trains(i_spiking_comp,:),2)
                    n_spiking_comp_accurate = n_spiking_comp_accurate + 1;
		    accurate(i_spiking_comp) = 1;
                end
            end
            if n_spikes_orig - n_spikes_comp == 1
                if spike_index(1) - compare_interval_sec*T.sampling_rate - best_match(i_spiking_comp,5) - 2*best_match(i_spiking_comp,6) < 0 || spike_index(end) + compare_interval_sec*T.sampling_rate + best_match(i_spiking_comp,5) + 2*best_match(i_spiking_comp,6) > size(spike_trains(i_spiking_comp,:),2)
                    n_spiking_comp_accurate = n_spiking_comp_accurate + 1;
		    accurate(i_spiking_comp) = 1;
                end
            end
        % every real spike is found in the sorted data, but sorted spike(s) is unaccounted for
        elseif size(find(curr_matches(:,4) == 1),1) == 1
            best_match(i_spiking_comp,:) = curr_matches(curr_matches(:,4) == 1,:);
            n_spikes_orig = round(best_match(i_spiking_comp,2)*best_match(i_spiking_comp,3));
            n_spikes_comp = best_match(i_spiking_comp,2);
            spike_index = find(spike_trains(i_spiking_comp,:));
            if n_spikes_comp - n_spikes_orig == 2
                if spike_index(1) - compare_interval_sec*T.sampling_rate - best_match(i_spiking_comp,5) - 2*best_match(i_spiking_comp,6) < 0 && spike_index(end) + compare_interval_sec*T.sampling_rate + best_match(i_spiking_comp,5) + 2*best_match(i_spiking_comp,6) > size(spike_trains(i_spiking_comp,:),2)
                    n_spiking_comp_accurate = n_spiking_comp_accurate + 1;
		    accurate(i_spiking_comp) = 1;
                end
            end
            if n_spikes_comp - n_spikes_orig == 1
                if spike_index(1) - compare_interval_sec*T.sampling_rate - best_match(i_spiking_comp,5) - 2*best_match(i_spiking_comp,6) < 0 || spike_index(end) + compare_interval_sec*T.sampling_rate + best_match(i_spiking_comp,5) + 2*best_match(i_spiking_comp,6) > size(spike_trains(i_spiking_comp,:),2)
                    n_spiking_comp_accurate = n_spiking_comp_accurate + 1;
		    accurate(i_spiking_comp) = 1;
                end
            end
        elseif size(find(curr_matches(:,3) == 1),1) > 1
            curr_matches_puta_exact = curr_matches(curr_matches(:,3) == 1,:);
            [val_max_real,i_max_real] = max(curr_matches_puta_exact(:,4));
            best_match(i_spiking_comp,:) = curr_matches_puta_exact(i_max_real,:);
        else
            [val_max_puta,i_max_puta] = max(curr_matches(:,3));
            best_match(i_spiking_comp,:) = curr_matches(i_max_puta,:);
        end
    end
    rate_false_pos_spikes(i_spiking_comp) = 1 - best_match(i_spiking_comp,3);
    rate_false_neg_spikes(i_spiking_comp) = 1 - best_match(i_spiking_comp,4);
    n_false_pos_spikes(i_spiking_comp) = (1 - best_match(i_spiking_comp,3))*best_match(i_spiking_comp,2);
    n_false_neg_spikes(i_spiking_comp) = best_match(i_spiking_comp,2)/best_match(i_spiking_comp,4) - best_match(i_spiking_comp,2);
end
disp(sprintf('%s %d matches found!',datestr(now),n_spiking_comp_accurate))