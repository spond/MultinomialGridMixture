_conditionals   = PATH_TO_CURRENT_BF + "data/conditional.txt";
_outputs        = PATH_TO_CURRENT_BF + "data/mcmc_samples.0";

fscanf (_outputs, "NMatrix,NMatrix", discard, weights);

_sample_count = Rows (weights);
cat_count  = Columns (weights);

prior_weights = {1, cat_count};

for (k = 0; k < _sample_count; k+=1) {
    prior_weights += weights[k][-1];
}

prior_weights = Transpose(prior_weights) * (1/_sample_count); 
prior_weights = prior_weights * (1/(+prior_weights));

conditionals = readMatrixFromText (_conditionals);
//fprintf ("data/dump.txt", CLEAR_FILE, conditionals);
//fscanf ("data/dump.txt", "NMatrix", conditionals);

_point_count = Rows (conditionals) ;
conditionals = conditionals [{{0,0}}][{{_point_count-1,cat_count-1}}];

// (i,j) == i-th site, j-th cat


normalization = conditionals * prior_weights;

posteriors = conditionals*({cat_count, cat_count}["(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)*prior_weights[_MATRIX_ELEMENT_ROW_]"]);
posteriors = ({_point_count, _point_count}["(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)*(1/normalization[_MATRIX_ELEMENT_ROW_])"])*posteriors;
 
/*
fprintf (stdout, Rows (posteriors), "x", Columns (posteriors), "\n"); 

for (k = 0; k < _point_count; k+=1) {
    fprintf (stdout, +(posteriors[k][-1]), ":", normalization[k], "\n");
}
*/

fprintf ("data/post.txt", CLEAR_FILE, Transpose(posteriors));

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