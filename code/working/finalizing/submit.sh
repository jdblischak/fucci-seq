inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(5; seq 10 10 100); do
 for ind in ${inds[@]}; do
#  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.ind.sbatch 15 $ngenes $ind
#  done
 done
done



inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
 for ind in ${inds[@]}; do
  sbatch job_run_methods.cyclical.ind.sbatch 15 $ind
done
