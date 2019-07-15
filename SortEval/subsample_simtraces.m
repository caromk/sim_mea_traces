n_electrode_rows = 321;
n_electrode_columns = 21;
dist_betw_electrodes = 1.25;

n_electrode = n_electrode_rows * n_electrode_columns;
            
x_coord_probe_edge = dist_betw_electrodes*(n_electrode_columns-1)/2;
            
% probe pad locations, centered on the rectangular area of the cylinder with z=0
loc_electrode = [reshape(repmat(-x_coord_probe_edge + [0:dist_betw_electrodes:dist_betw_electrodes*(n_electrode_columns-1)]',1,n_electrode_rows)',1,[]);...
    repmat(-dist_betw_electrodes*(n_electrode_rows-1)/2 + [0:dist_betw_electrodes:dist_betw_electrodes*(n_electrode_rows-1)],1,n_electrode_columns);...
    zeros(1,n_electrode)]';

loc_index = ismember(CT.loc_electrode,loc_electrode,'rows');