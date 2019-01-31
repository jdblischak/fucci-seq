#---- identify top cylcial in each training set and compute prediction error in each fold
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
 for ind in ${inds[@]}; do
  sbatch job_run_methods.cyclical.ind.sbatch 15 $ind
done

inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 5 5; seq 10 10 500); do
 for ind in ${inds[@]}; do
  sbatch job_run_methods.train.ind.sbatch 15 $ngenes $ind
 done
done


# --- use the top cylical across the entire dataset to compute prediction error in each fold
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 5 5; seq 10 10 50); do
 for ind in ${inds[@]}; do
  sbatch job_run_methods.train.ind.overall.sbatch 15 $ngenes $ind
 done
done

inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 210 10 500); do
 for ind in ${inds[@]}; do
  sbatch job_run_methods.train.ind.overall.sbatch 15 $ngenes $ind
 done
done



# # make cyclone pairs gene using Leng data
# sbatch job_run_methods.cyclone.sbatch
#
#
#
# # Do mixed indivdiual
# for ngenes in $(seq 5 5; seq 10 10 20); do
#  for fold in $(seq 1 6); do
#   sbatch job_run_methods.train.overall.sbatch 15 $ngenes $fold
#  done
# done
#
# for ngenes in $(seq 30 10 500); do
#  for fold in $(seq 1 6); do
#   sbatch job_run_methods.train.overall.sbatch 15 $ngenes $fold
#  done
done

