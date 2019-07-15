function CT = subsetT_by_params(CT,n_row,n_col,pitch,varargin)

% parse inputs
p = inputParser;
default_diamond = 0;
addOptional(p,'diamond',default_diamond,@isnumeric);
p.parse(varargin{:});

x_locs = -(n_col-1)/2*pitch:pitch:(n_col-1)/2*pitch;
y_locs = -(n_row-1)/2*pitch:pitch:(n_row-1)/2*pitch;

ind_use = find(ismember(CT.loc_electrode(:,1),x_locs) & ismember(CT.loc_electrode(:,2),y_locs));

if p.Results.diamond
    n_row2 = n_row-1;
    n_col2 = n_col-1;
    x_locs2 = -(n_col2-1)/2*pitch:pitch:(n_col2-1)/2*pitch;
    y_locs2 = -(n_row2-1)/2*pitch:pitch:(n_row2-1)/2*pitch;
    ind_use2 = find(ismember(CT.loc_electrode(:,1),x_locs2) & ismember(CT.loc_electrode(:,2),y_locs2));
    ind_use = [ind_use;ind_use2];
end

[val_unique ind_unique] = unique(CT.loc_electrode(ind_use,:),'rows');
CT = subsetT(CT,ind_use(ind_unique));
