function CT = subset(CT,row,col,pitch)

x_locs = -(col-1)/2*pitch:pitch:(col-1)/2*pitch;
y_locs = -(row-1)/2*pitch:pitch:(row-1)/2*pitch;

ind_use = find(ismember(CT.loc_electrode(:,1),x_locs) & ismember(CT.loc_electrode(:,2),y_locs));

[val_unique ind_unique] = unique(CT.loc_electrode(ind_use,:),'rows');
CT = subsetT(CT,ind_use(ind_unique));