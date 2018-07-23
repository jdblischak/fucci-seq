# Time basd on fucci --------------------------------------------------------------

# Get cyclical PVE

# DONE
for fold in $(seq 1 5); do
  sbatch job_run_methods.cyclical.sbatch 20 $fold
done

# DONE
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.cyclical.ind.sbatch 15 $fold $ind
  done
done


# Prediction error

# done
for ngenes in $(seq 10; seq 20 10 800); do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.sbatch 15 $ngenes $fold
  done
done

# submitted
for ngenes in $(seq 810 10 1000); do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.sbatch 15 $ngenes $fold
  done
done


# DONE
for ngenes in $(seq 5; seq 10 10 200); do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.permute.oncyclical.sbatch 15 $ngenes $fold
  done
done


# done
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 9; seq 10 10 300); do
 for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.ind.sbatch 15 $ngenes $fold $ind
  done
 done
done

# done?
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 310 10 500); do
 for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.ind.sbatch 15 $ngenes $fold $ind
  done
 done
done

# submitted (?)
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 510 10 1000); do
 for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.ind.sbatch 15 $ngenes $fold $ind
  done
 done
done


inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 410 10 490); do
 for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.ind.sbatch 15 $ngenes $fold $ind
  done
 done
done

#TBD
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 680 10 800); do
 for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.ind.sbatch 15 $ngenes $fold $ind
  done
 done
done


# TBD
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 5 5 ; seq(10 20 200); do
 for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.permute.oncyclical.ind.sbatch 15 $ngenes $fold $ind
  done
 done
done






# Time based on fucci and DAPI ---------------------------------------------------------

# Get cyclical PVE

# done
for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.cyclical.sbatch 15 $fold
done

# done
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.cyclical.ind.sbatch 15 $fold $ind
  done
done


# Get prediction error

# done
for ngenes in $(seq 5 5; seq 10 10 200); do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.train.sbatch 15 $ngenes $fold
  done
done

# done
for ngenes in $(seq 5 5; seq 10 20 200); do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.train.permute.oncyclical.sbatch 15 $ngenes $fold
  done
done

# TBD
for ngenes in $(seq 210 20 400); do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.train.permute.oncyclical.sbatch 15 $ngenes $fold
  done
done

# TBD
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 9); do
for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.train.ind.sbatch 15 $ngenes $fold $ind
  done
done
done



#TBD
for ngenes in $(seq 9; seq 10 10 200); do
for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.train.permute.oncyclical.ind.sbatch 15 $ngenes $fold $ind
  done
done
done


# triple time based on algebraic methods ------------------------------------------

# Get cyclical PVE

# submitted
for fold in $(seq 1 5);
do
  sbatch job_run_methods.triple.cyclical.alge.sbatch 15 $fold
done

# submitted
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.cyclical.ind.alge.sbatch 15 $fold $ind
  done
done


# Prediction error

# done
for ngenes in $(seq  2 9; seq 10 20 200); do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.train.alge.sbatch 15 $ngenes $fold
  done
done



# done
inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in $(seq 2 9; seq 10 10 200); do
for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.train.ind.alge.sbatch 15 $ngenes $fold $ind
  done
done
done

# # time based on fucci using algebraic methods --------------------------------------------
#
# # done
# for fold in $(seq 1 5);
# do
#   sbatch job_run_methods.cyclical.alge.sbatch 20 $fold
# done
#
# # for ngenes in 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200; do
# #   for fold in $(seq 1 5); do
# #   sbatch job_run_methods.train.sbatch 15 $ngenes $fold
# #   done
# # done
#
#
# # # compute prediction error of top X genes for the random cell times
# # # where the top genes are selected from the observed fucci data
# # for ngenes in 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200; do
# #   for fold in $(seq 1 5); do
# #   sbatch job_run_methods.train.permute.oncyclical.sbatch 15 $ngenes $fold
# #   done
# # done
#
# # # submitted
# # inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
# # for ind in ${inds[@]}; do
# #   for fold in $(seq 1 5); do
# #   sbatch job_run_methods.cyclical.ind.alge.sbatch 15 $fold $ind
# #   done
# # done
#
#
# # # TBD
# # inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
# # for ngenes in 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200; do
# # # for ind in ${inds[@]}; do
# #   for fold in $(seq 1 5); do
# #   sbatch job_run_methods.train.ind.sbatch 15 $ngenes $fold $ind
# #   done
# # # done
# # done
# #
# # # TBD
# # for ngenes in 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200; do
# # # for ind in ${inds[@]}; do
# #   for fold in $(seq 1 5); do
# #   sbatch job_run_methods.train.permute.oncyclical.ind.sbatch 15 $ngenes $fold $ind
# #   done
# # # done
# # done
#





