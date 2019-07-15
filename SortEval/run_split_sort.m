% wrapper for compiled matlab code to run spatially split spike sorting on
% sim data
function run_split_sort(filename,split_row,split_col)

if isa(split_row,'char')
    split_row = str2double(split_row);
end

if isa(split_col,'char')
    split_col = str2double(split_col);
end 

load(filename);

[E comps all_thres_spikes spike_comps wts] = SortEval_cruns(CT,'SplitRow',split_row,'SplitCol',split_col);

neurons_steps = 10:10:200;
return_ratio = zeros(1,length(neurons_steps));
for i = 1:length(neurons_steps)
    return_ratio(i) = sum(E.min_dist_from_probe_center < neurons_steps(i))/sum(dist_from_probe_center < neurons_steps(i));
end

save('-v7.3',sprintf('%s_split.mat',filename),'E','comps','all_thres_spikes','spike_comps','wts','neurons_steps','return_ratio','split_row','split_col');
