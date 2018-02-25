%MS_batch_control

job = batch('MS_discovery_batch', 'Pool', 7);
wait(job)
load(job, 'mat_all')

save(['/dartfs-hpc/rc/lab/M/MeerM/EC/temp/' 'MS_mat2.mat'], 'mat_all', '-v7.3')

