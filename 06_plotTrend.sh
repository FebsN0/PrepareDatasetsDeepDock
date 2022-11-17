#!/bin/bash
#SBATCH --account=def-jtus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --gres=gpu:1
#SBATCH --mem=5G               # memory per node
#SBATCH --job-name=plotTrend
#SBATCH --time=1:0:0

# $1 : it_start
# $2 : it_end
if [[ ! -n $1 || ! -n $2 ]]; then echo -e "\n\tstart and/or end iteration missing\n"; exit 1; fi;
# $3: file logs.txt
if [[ ! -n $3 ]]; then echo "logs.txt not loaded!"; exit 1; fi;

source ~/envDeepDock/bin/activate

file_path=`sed -n '1p' $3`
protein=`sed -n '2p' $3`    # name of project folder
sz_test=`sed -n '6p' $3`
if [[ -d resultTrend_"$1"_"$2" ]]; then rm -r resultTrend_"$1"_"$2"; fi;
python $DEEPDOCKNA/utilities/plot_progress.py -path_protein $file_path/$protein -sz_test $sz_test -it_start $1 -it_end $2 -fo resultTrend_"$1"_"$2"
