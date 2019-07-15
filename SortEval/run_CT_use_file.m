function run_CT_use_file(filename)

load(sprintf('%s%s.mat',filename));
[E_use comps_use] = SortEval_cruns(CT_use);
save('-v7.3',sprintf('%s_comps.mat',filename),'E_use','comp_use');