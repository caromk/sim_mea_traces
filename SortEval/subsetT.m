function Tsub = subsetT(T,sub_ind)

Tsub = T;
Tsub.traces = Tsub.traces(sub_ind,:);
Tsub.dist_from_soma = Tsub.dist_from_soma(sub_ind,:);
Tsub.loc_electrode = Tsub.loc_electrode(sub_ind,:);
Tsub.n_electrode = length(sub_ind);
