using DelimitedFiles

batch_data_dir = "../../data/phase_diag_N300/g_rel_values/"
mkpath(batch_data_dir)

#number of batches 
nr_batch = 3

nr_vals = 80
g_rel_values = range(0.1,10,length=nr_vals) 
batch_size = nr_vals ÷ nr_batch #won't be of equal size (last one is different)
batch_start_idxs = 1:batch_size:nr_vals

for i in 1:length(batch_start_idxs)-1
    if i !== nr_batch 
        batch_vals = g_rel_values[batch_start_idxs[i]:batch_start_idxs[i+1]-1]
    else 
        batch_vals = g_rel_values[batch_start_idxs[i]:end]
    end
    
    writedlm(batch_data_dir * "g_rel_vals_batch$(i).txt" ,batch_vals)
    
end
