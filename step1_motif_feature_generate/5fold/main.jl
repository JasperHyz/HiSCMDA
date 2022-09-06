#=
run:
- Julia version:
- Author: Allen
- Date: 2018-08-17
=#
inputfile= "motif_step1_SEED_0_FOLD_0"
outputfile= "motif_step2_SEED_0_FOLD_0"
include("util.jl");
para = algPara(0.8, 5, 100, 0.4, outputfile);
P = read_tensor(inputfile);
tensor_srw(P, para);
