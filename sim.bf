MESSAGE_LOGGING = 0;
LoadFunctionLibrary ("ProbabilityDistributions");

sampleFromThisDistro = Transpose ({{0,1,2,3}
                        {0.98,0.015,0.025,0.025}});

for (site = 0; site < 500; site += 1) {
    sims = Random(sampleFromThisDistro, {"PDF":"Multinomial","ARG0":sampleIndependentNormal (10000,1000)});
    fprintf (stdout, Join (" ", sims[-1][1]), "\n");
}

lfunction sampleIndependentNormal (mu, sigma) {
    return sampleFromNormal () * sigma + mu;
}  