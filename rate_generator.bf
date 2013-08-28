majority_freqs = {{0.999, 0.995, 0.99, 0.95, 0.9, 0.8, 0.75, 0.5, 0.33}};
mj_fc = Columns (majority_freqs);
minority_freqs = {{0.01/3, 0.005/3, 0.025, 0.005, 0.01, 0.05/3, 0.1/3, 0.05, 0.10, 0.25}};
mn_fc = Columns (minority_freqs);


for (majority_residue = 0; majority_residue < 4; majority_residue += 1) {
    for (mj_freq = 0; mj_freq < mj_fc; mj_freq += 1) {
        freqs = {1,4};
        freqs[majority_residue] = majority_freqs[mj_freq];
        minority_configuration = {3,1};
        minority_map = {3,1};
        for (k = 0; k < majority_residue; k+=1) {
            minority_map[k] = k;
        }
        for (k = majority_residue+1; k < 4; k+=1) {
            minority_map[k-1] = k;
        }
        
        for (r1 = 0; r1 <= mn_fc; r1 += 1) {
            minority_configuration[0] = r1; 
            for (r2 = 0; r2 <= mn_fc; r2 += 1) {
                minority_configuration[1] = r2;
                for (r3 = 0; r3 <= mn_fc; r3 += 1) {
                    minority_configuration[2] = r3;
                    // validate the configuration 
                    for (k = 0; k < 3; k += 1) {
                        if (minority_configuration[k]) {
                            freqs[minority_map[k]] = minority_freqs[minority_configuration[k]-1];
                        }
                    }
                    total_set = (+minority_configuration["_MATRIX_ELEMENT_VALUE_>0"]);
                    if (total_set < 3) {
                        if (+freqs < 1) {
                            left_over = (1-(+freqs))/(3-total_set);
                            freqs = freqs["_MATRIX_ELEMENT_VALUE_*(_MATRIX_ELEMENT_VALUE_>0)+(_MATRIX_ELEMENT_VALUE_==0)*left_over"];
                            fprintf (stdout, Join (" ", freqs[{{0,0}}][{{0,2}}]), "\n");
                        }
                    } else {
                        if (+freqs == 1) {
                            fprintf (stdout, Join (" ", freqs[{{0,0}}][{{0,2}}]), "\n");                        
                        }
                    }
                }
            }
        }
    }
}