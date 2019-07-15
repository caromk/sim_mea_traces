function [E,comps,varargout] = SortEval_cruns2(T,varargin)

default_alt_length = size(T.traces,2);
default_comps = [];
default_verbose = 0;
default_max_peak_variation = 0.4;
default_min_separation = 0.3;
default_compare_interval_sec = 1.5e-3;
default_use_binary = 0;
default_split_row = 0;
default_split_col = 0;
default_split_time = 0;

% get options 
opt = Opt(varargin);

% can run with a shorter length
alt_length = opt.get('AltLength',default_alt_length);
% init optional inputs / defaults
% verbose prints out details and displays figures
verbose = opt.get('Verbose',default_verbose);
% ratio for allowed variation in peak height, in ratio of peak height
max_peak_variation = opt.get('MaxPeakVariation',default_max_peak_variation);
% minimum separation distance (in max_peak scale) between nonzero intervals
% for the last nonzero interval to be the robust peaks
min_separation = opt.get('MinSeparation',default_min_separation);
% compare interval is the window used to compare spikes and spike trains to each other to see if they're the same 
compare_interval_sec = opt.get('CompareIntervalSec',default_compare_interval_sec);
% use binary file to run ICA (otherwise, run matlab version)
use_binary = opt.get('UseBinary',default_use_binary);
% subset into columns, number of columns per subset
split_col = opt.get('SplitCol',default_split_col);
% subset into rows, number of rows per subset
split_row = opt.get('SplitRow',default_split_row);
% subset in time, seconds per subset
split_time = opt.get('SplitTime',default_split_time);
% comps = components returned from ICA, added as an option for cases in
% which the algorithm needs to be run multiple times, makes it faster and
% more consistent
comps = opt.get('Comps',default_comps);
if isempty(comps)
    fprintf('Running comps.\n')
end

% prep row / col split if using
ind_space_splits = {};
if split_col ~= 0 || split_row ~=0
    uniq_y = sort(unique(T.loc_electrode(:,2)));
    n_rows = length(uniq_y);
    if split_row == 0
        split_row = n_rows;
    end
    uniq_x = sort(unique(T.loc_electrode(:,1)));
    n_cols = length(uniq_x);
    if split_col == 0
        split_col = n_cols;
    end
    tile_rows = ceil(2*n_rows/split_row-1);
    tile_cols = ceil(2*n_cols/split_col-1);
    tile_row_pitch = (n_rows - split_row)/(tile_rows-1);
    tile_col_pitch = (n_cols - split_col)/(tile_cols-1);
    i_tile = 0;
    for i_tile_row = 1:tile_rows
        for i_tile_col = 1:tile_cols
            i_tile = i_tile + 1;            
            if tile_cols == 1
                x_values = uniq_x;
            else
                x_start = ceil((i_tile_col-1)*tile_col_pitch+1);
                if x_start > 1 && x_start + split_col - 1 > n_cols
                    x_start = n_cols - split_col + 1;
                end
                x_values = uniq_x(x_start:(x_start+split_col-1));
            end
            if tile_rows == 1
                y_values = uniq_y;
            else
                y_start = ceil((i_tile_row-1)*tile_row_pitch+1);
                if y_start > 1 && y_start + split_row - 1 > n_rows
                    y_start = n_rows - split_row + 1;
                end
                y_values = uniq_y(y_start:(y_start+split_row-1));
            end
            ind_space_splits{i_tile} = find(ismember(T.loc_electrode(:,1),x_values) & ismember(T.loc_electrode(:,2),y_values));   
        end
    end
end

tic
[E.spike_trains,spike_comps,comps,E.n_dups,E.best_spike_index,E.ind_spiking_comps,E.separation,E.peak_variation,wts,dups,E.separation_noise,E.peak_variation_noise,E.dist_from_noise] = RobustSpikeSort2(T.traces(:,1:alt_length),T.sampling_rate,'Verbose',verbose,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'Comps',comps,'IndSpaceSplits',ind_space_splits,'UseBinary',use_binary);
E.run_time = toc;
E.max_peak_variation = max_peak_variation;
E.min_separation = min_separation;

E.n_spiking_comp = size(E.ind_spiking_comps,2);
E.matches = [];
E.best_match = zeros(E.n_spiking_comp,6);
E.n_spiking_comp_accurate = 0;
E.rate_false_pos_spikes = zeros(E.n_spiking_comp,1);
E.rate_false_neg_spikes = zeros(E.n_spiking_comp,1);
E.n_false_pos_spikes = zeros(E.n_spiking_comp,1);
E.n_false_neg_spikes = zeros(E.n_spiking_comp,1);
E.accurate = zeros(E.n_spiking_comp,1);

disp(sprintf('%s %d spiking components',datestr(now),E.n_spiking_comp))
[E.matches,E.best_match,E.n_spiking_comp_accurate,E.rate_false_pos_spikes,E.rate_false_neg_spikes,E.n_false_pos_spikes,E.n_false_neg_spikes,E.accurate] = EvalSorted(T,E.spike_trains,compare_interval_sec);

%E.n_neurons_50microV = sum(sum(T.calc_snr_multiplier(T.dist_from_soma)>50)>=1);
E.min_dist_from_pad = zeros(E.n_spiking_comp,1);
E.min_dist_from_pad(E.best_match(:,1)>=1) = min(T.dist_from_soma(:,E.best_match(E.best_match(:,1)>=1,1)))';

y_max = max(T.loc_electrode(:,2));
E.min_dist_from_probe_center_line = zeros(E.n_spiking_comp,1);
goodmatches = sum(E.best_match(:,1)>=1);
E.min_dist_from_probe_center_line(E.best_match(:,1)>=1) = sqrt(T.loc_neuron(E.best_match(E.best_match(:,1)>=1,1),1).^2 + max(zeros(goodmatches,1),abs(T.loc_neuron(E.best_match(E.best_match(:,1)>=1,1),2))-y_max).^2 + T.loc_neuron(E.best_match(E.best_match(:,1)>=1,1),3).^2);
E.min_dist_from_probe_center = sqrt(sum(T.loc_neuron(E.best_match(E.best_match(:,1)>=1,1),:).^2,2));

varargout{1} = spike_comps;
%varargout{2} = wts;
