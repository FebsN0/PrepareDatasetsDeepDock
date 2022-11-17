#!/bin/bash
#SBATCH --account=def-jtus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --gres=gpu:1
#SBATCH --mem=10G               # memory per node
#SBATCH --job-name=phase_5
#SBATCH --time=2:0:0

# $1 : n_it
if [[ ! -n $1 ]]; then echo -e "\n\tnumber of current iteration missing\n"; exit 1; fi;
# $2: file logs.txt
if [[ ! -n $2 ]]; then echo "logs.txt not loaded!"; exit 1; fi;

source ~/envDeepDock/bin/activate

file_path=`sed -n '1p' $2`
protein=`sed -n '2p' $2`    # name of project folder

morgan_directory=`sed -n '3p' $2`
num_molec=`sed -n '6p' $2`


echo "Starting Evaluation"
python -u $DEEPDOCKNA/phase_2-3/hyperparameter_result_evaluation.py -n_it $1 -d_path $file_path/$protein -mdd $morgan_directory -n_mol $num_molec
echo "Creating simple_job_predictions"

python $DEEPDOCKNA/phase_2-3/simple_job_predictions.py -protein $protein -file_path $file_path -n_it $1 -mdd $morgan_directory

cd $file_path/$protein/iteration_$1/simple_job_predictions/
echo "running simple_jobs"
for f in *;do sbatch $f;done
