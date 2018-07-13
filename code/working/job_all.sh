for fold in $(seq 1 5);
do
  sbatch job_run_methods.cyclical.sbatch 20 $fold
done

for ngenes in 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.sbatch 15 $ngenes $fold
  done
done

for ngenes in 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.sbatch 15 $ngenes $fold
  done
done

# compute PVE under random cell times
for fold in $(seq 1 5); do
  sbatch job_run_methods.cyclical.permute.sbatch 15 $fold
done

# compute prediction error of top X genes selected from random cell times
for ngenes in 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.permute.sbatch 15 $ngenes $fold
  done
done

# for ngenes in 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400; do
#   for fold in $(seq 1 5); do
#   sbatch job_run_methods.train.permute.sbatch 15 $ngenes $fold
#   done
# done

# compute prediction error of top X genes for the random cell times
# where the top genes are selected from the observed fucci data
for ngenes in 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.permute.oncyclical.sbatch 15 $ngenes $fold
  done
done


for ngenes in 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.permute.onrandom.sbatch 15 $ngenes $fold
  done
done

### validation sample

sbatch job_run_methods.eval.seuratgenes.sbatch 15

for nlist in $(seq 1 41); do
#nlist=1
  sbatch job_run_methods.eval.sbatch 15 $nlist
done


for nlist in $(seq 1 41); do
#nlist=1
  sbatch job_run_methods.eval.permute.oncyclical.sbatch 15 $nlist
done



############# validation on ranodm individual

inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")

for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.cyclical.ind.sbatch 15 $fold $ind
  done
done

#inds=("NA19098" "NA18511" "NA18870" "NA19101" "NA18855" "NA19160")
for ngenes in 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200; do
ind="NA18511"
# for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.train.ind.sbatch 15 $ngenes $fold $ind
  done
# done
done



for ind in ${inds[@]}; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.cyclical.ind.permute.sbatch 15 $fold $ind
# #  sbatch job_run_methods.train.ind.sbatch 15 $ngens $fold $ind
  done
done





############# triple time
for fold in $(seq 1 5);
do
  sbatch job_run_methods.triple.cyclical.sbatch 15 $fold
done

for ngenes in 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200; do
  for fold in $(seq 1 5); do
  sbatch job_run_methods.triple.train.sbatch 15 $ngenes $fold
  done
done
