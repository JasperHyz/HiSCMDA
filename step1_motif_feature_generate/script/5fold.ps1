conda activate python38
cd ../5fold

# set main.jl origin
# origin main.jl
# inputfile= "motif_step1_SEED_0_FOLD_0"
# outputfile= "motif_step2_SEED_0_FOLD_0"
$SEED=0

for($FOLD=1;$FOLD -le 5;$FOLD++)
{
    python triangle2miRNA_MotifGenerator.py -s $SEED -f $FOLD
    python triangle2miRNA_TensorGenerator.py -s $SEED -f $FOLD
    $file = 'main.jl'
    $motif = 'triangle2miRNA'
    $find1 = "motif_step1_SEED_0_FOLD_0"
    $replace1 = "{0}_step1_SEED_{1}_FOLD_{2}" -f $motif,$SEED,$FOLD
    $find2 = "motif_step2_SEED_0_FOLD_0"
    $replace2 = "{0}_step2_SEED_{1}_FOLD_{2}" -f $motif,$SEED,$FOLD
    (Get-Content $file).replace($find1, $replace1) | Set-Content $file
    (Get-Content $file).replace($find2, $replace2) | Set-Content $file
    # running julia 
    julia main.jl
    # back to origin
    (Get-Content $file).replace($replace1, $find1) | Set-Content $file
    (Get-Content $file).replace($replace2, $find2) | Set-Content $file
    python post_process.py -s $SEED -f $FOLD -m $motif


    python triangle2disease_MotifGenerator.py -s $SEED -f $FOLD
    python triangle2disease_TensorGenerator.py -s $SEED -f $FOLD
    $motif = 'triangle2disease'
    $replace1 = "{0}_step1_SEED_{1}_FOLD_{2}" -f $motif,$SEED,$FOLD
    $replace2 = "{0}_step2_SEED_{1}_FOLD_{2}" -f $motif,$SEED,$FOLD
    (Get-Content $file).replace($find1, $replace1) | Set-Content $file
    (Get-Content $file).replace($find2, $replace2) | Set-Content $file
    # running julia 
    julia main.jl
    # back to origin
    (Get-Content $file).replace($replace1, $find1) | Set-Content $file
    (Get-Content $file).replace($replace2, $find2) | Set-Content $file
    python post_process.py -s $SEED -f $FOLD -m $motif


    python triangle_MotifGenerator.py -s $SEED -f $FOLD
    python triangle_TensorGenerator.py -s $SEED -f $FOLD
    $motif = 'triangle'
    $replace1 = "{0}_step1_SEED_{1}_FOLD_{2}" -f $motif,$SEED,$FOLD
    $replace2 = "{0}_step2_SEED_{1}_FOLD_{2}" -f $motif,$SEED,$FOLD
    (Get-Content $file).replace($find1, $replace1) | Set-Content $file
    (Get-Content $file).replace($find2, $replace2) | Set-Content $file
    # running julia 
    julia main.jl
    # back to origin
    (Get-Content $file).replace($replace1, $find1) | Set-Content $file
    (Get-Content $file).replace($replace2, $find2) | Set-Content $file
    python post_process.py -s $SEED -f $FOLD -m $motif


    python mdm_MotifGenerator.py -s $SEED -f $FOLD
    python mdm_TensorGenerator.py -s $SEED -f $FOLD
    $motif = 'mdm'
    $replace1 = "{0}_step1_SEED_{1}_FOLD_{2}" -f $motif,$SEED,$FOLD
    $replace2 = "{0}_step2_SEED_{1}_FOLD_{2}" -f $motif,$SEED,$FOLD
    (Get-Content $file).replace($find1, $replace1) | Set-Content $file
    (Get-Content $file).replace($find2, $replace2) | Set-Content $file
    # running julia 
    julia main.jl
    # back to origin
    (Get-Content $file).replace($replace1, $find1) | Set-Content $file
    (Get-Content $file).replace($replace2, $find2) | Set-Content $file
    python post_process.py -s $SEED -f $FOLD -m $motif

    
    python mddd_MotifGenerator.py -s $SEED -f $FOLD
    python mddd_TensorGenerator.py -s $SEED -f $FOLD
    $motif = 'mddd'
    $replace1 = "{0}_step1_SEED_{1}_FOLD_{2}" -f $motif,$SEED,$FOLD
    $replace2 = "{0}_step2_SEED_{1}_FOLD_{2}" -f $motif,$SEED,$FOLD
    (Get-Content $file).replace($find1, $replace1) | Set-Content $file
    (Get-Content $file).replace($find2, $replace2) | Set-Content $file
    # running julia 
    julia main.jl
    # back to origin
    (Get-Content $file).replace($replace1, $find1) | Set-Content $file
    (Get-Content $file).replace($replace2, $find2) | Set-Content $file
    python post_process.py -s $SEED -f $FOLD -m $motif
}