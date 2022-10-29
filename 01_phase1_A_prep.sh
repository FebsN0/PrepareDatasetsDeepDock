#!/bin/bash

#SBATCH --mem=10G
#SBATCH --account=def-jtus
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00

# $1: file logs.txt
# $2: current iteration

#activate virtual environment where there all installed packages
source ~/envDeepDock/bin/activate

#read -p 'current iteration? ' currIt
#read -p 'number of CPUs (cores)? ' t_cpu

t_cpu=$SLURM_CPUS_PER_TASK
currIt=$2

# main path where there is anything
file_path=`sed -n '1p' $1`
# name project: ie name of directory where save anything
protein=`sed -n '2p' $1`
mkdir $protein
# name dir where there is morgan FP data
morgan_directory=`sed -n '3p' $1`
# name dir where there is SMILE data
smile_directory=`sed -n '4p' $1`
# size of initial dataset (suggested 1milion for train, test, validation sets)
n_mol=`sed -n '5p' $1`

trainSize=$(($n_mol/3))
validTestSize=$trainSize


# current iteration
pr_it=$(($currIt-1))

if [ $currIt == 1 ]
then
	pred_directory=$morgan_directory
else
	pred_directory=$file_path/$protein/iteration_$pr_it/morgan_1024_predictions
fi


# $DEEPDOCKNA is where the DeepDock-NonAutomated directry is located. Exported in .bashrc

#generate in $pred_directory (if it=1, it's morganDirectory, otherwise inside protein dir) 2 similar files: Mol_ct_file.csv and Mol_ct_file_update.csv
python $DEEPDOCKNA/phase_1/molecular_file_count_updated.py -pt $protein -it $currIt -cdd $pred_directory -t_pos $t_cpu -t_samp $n_mol
echo -e 'END MOLECULAR COUNT\n\n\n'

####ADD TRAINING and VALID/TEST SIZE
#generate training, test, validition datasets. Each line is ZINCID only. Based on number of given trainSize validTestSize
python $DEEPDOCKNA/phase_1/sampling.py -pt $protein -fp $file_path -it $currIt -dd $pred_directory -t_pos $t_cpu -tr_sz $trainSize -vl_sz $validTestSize
echo -e 'END RANDOM SAMPLING\n\n\n'

#remove overlap sampled sets, checking training/validation/test_labels*
python $DEEPDOCKNA/phase_1/sanity_check.py -pt $protein -fp $file_path -it $currIt
echo -e 'END ty CHECK\n\n\n'

#extract morganFP and SMILE compound from the directory where there are all bilion compounds. Save in protein/iteration$i/{morgan}{smile}
python $DEEPDOCKNA/phase_1/Extracting_morgan.py -pt $protein -fp $file_path -it $currIt -md $morgan_directory -t_pos $t_cpu
echo -e 'END EXTRACTION MORGAN\n\n\n'

python $DEEPDOCKNA/phase_1/Extracting_smiles.py -pt $protein -fp $file_path -it $currIt -smd $smile_directory -t_pos $t_cpu
echo -e 'END EXTRACTION SMILE\n\n\n'


