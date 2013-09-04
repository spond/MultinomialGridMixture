_posteriors     = PATH_TO_CURRENT_BF + "data/post.txt";
_gridInfo       = PATH_TO_CURRENT_BF + "data/rates.txt";

threshold   = 0.1;
post_thresh = 0.999;

_rates = readMatrixFromText (_gridInfo);
rate_count = Rows (_rates);
passThresh = {rate_count, 4};

for (k = 0; k < rate_count; k += 1) {
    sum = 0;
    for (n = 0; n < 3; n+=1) {
        passThresh[k][n] = _rates[k][n] > threshold;
        sum += _rates[k][n];
    }
    passThresh[k][n] = (1-sum) >= threshold;
}


fscanf (_posteriors, "NMatrix", _posteriors);

site_count = Columns (_posteriors);
filtered_posteriors = {site_count, 4};

haz = Transpose (_posteriors) * (passThresh);
chars = "ACGT";

/*for (k = 0; k < rate_count; k += 1) {
    if (_posteriors[k][2717] > 0.001) {
        fprintf (stdout, _rates[k][2], " ", 1 - (+_rates[k][-1]), ":", _posteriors[k][2717], "\n");
    }
}*/


for (s = 0; s < site_count; s+=1) {
   fprintf (stdout, s+1, " ");
   for (n = 0; n < 4; n+=1) {
    if (haz[s][n] >= post_thresh) {
        fprintf (stdout, chars[n]);
    }
   }
   fprintf (stdout, "\n");
}

//------------------------------------------------------------------------------------------------//

lfunction readMatrixFromText (filename) {
    fscanf (filename, REWIND, "Lines", file_contents);
    column_counter = 1 + Rows (file_contents[0]||"\\ ") $ 2;
    row_counter = Columns (file_contents);
    matrix = {row_counter, column_counter};
    for (r = 0; r < row_counter; r += 1) {
        row = file_contents[r];
        sscanf (row, REWIND, "Number", v);
        matrix[r][0] = v;
        for (c = 1; c < column_counter; c += 1) {
            sscanf (row, "Number", v);
            /*if (c == column_counter - 1 && c > 3 && r > 0) {
                fprintf (stdout, column_counter, "\n", v, "\n", row, "\n");
                assert (0);
            }*/
            matrix[r][c] = v;
        }
    }
    file_contents = 0;
    return matrix;
}