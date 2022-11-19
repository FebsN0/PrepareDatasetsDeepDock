#!/bin/bash

# $1: file logs.txt
if [[ $1 == logs.txt ]]
then
        echo "logs.txt missing"
        exit 1
fi

file_path=`sed -n '1p' $1`
morgan_directory=`sed -n '3p' $1`

# name project: ie name of directory where save anything
protein=`sed -n '2p' $1`
read -p "name project is $protein, change or not? [0: change | other/none: ok] " ans
case $ans in
        0) read -p 'nameDirectoryOfProject: ' protein;;
        *) ;;
esac
if [ ! -d $file_path/$protein ]; then echo "$protein directory doesn't exist!"; fi

#current iteration
read -p "current iteration: " currIt
if [ ! -d $file_path/$protein/iteration_$currIt ]; then echo "Iteration $currIt dir inside $protein dir doesnt exist"; exit 1; fi;

cd $file_path/$protein/iteration_$currIt


#select the binding site found previously
echo -e "\nBINDING SITES AVAILABLE:\n 1 (data saved on graham)\n 2 (data saved on narval)\n 3 (data saved on cedar)\n 4 (data saved on beluga)"
read -p 'choose the binding site? ' siteSel
case $siteSel in
        1|2|3|4);;
        *) echo "wrong number. Select only 1,2,3 or 4"; exit 1;;
esac
echo -e ""
sets=(train test valid)

if [[ $currIt -eq 1 ]]
then
	echo -e "DATASETs:\n 0 (TRAINING SET - train Dir)\n 1 (TEST SET - test Dir)\n 2 (VALIDATION SET - valid Dir)"
	read -p 'which dataset to run? ' setSel
	case $setSel in
        	0|1|2);;
        	*) echo "wrong answer. Select only 0,1 or 2"; exit 1;;
	esac
else
	setSel=0
	cp ../iteration_1/testing_labels.txt ../iteration_1/validation_labels.txt .
fi

#parallelize some jobs, especially during the restart confSearch, starting/resuming docking jobs and postprocessing jobs
echo ""
read -p 'how many cpus use to parallelize? ' t_cpu


#the python script crash with more than 512 blocks from whuch extract the data.
cd docking/site_$siteSel/${sets[$setSel]}
if [[ -n $(ls) ]]
then
        lastBlockProc=$(ls | sort -n -t 'k' -k 2 | tail -n 1)
	lastBlock="${lastBlockProc##*k}"          #the format is blockNUMBER
	echo -e "number of processed blocks:\t$lastBlock"
else
        lastBlockProc=0
        lastBlock=0
        echo -e "\nThere are no processed blocks\nExiting because not possible to check and process. Run first the code 02\n"
        exit 0
fi
cd ..

if [[ ! -d dock_result_blocksChunked_${sets[$setSel]} ]]
then
	read -p "how many blocks put in macro chunks you want create? (having more than 400-500 blocks create problems) " numMacroC
	if [[ $numMacroC -lt 1 && $numMacroC -gt $lastRun ]]
	then
		echo "wrong number"
		exit 1
	fi
fi

#if the operation is not done yet. check if number of blocks are over 500, if yes,create macroblocks
if [[ $lastBlock -gt $numMacroC && ! -d dock_result_blocksChunked_${sets[$setSel]} ]]
then
	nMacroBlocks=$(($lastBlock/$numMacroC))	# number of macroBlock containing 500 blocks
	nRestBlocks=$(($lastBlock%$numMacroC))	# number of blocks rest
	echo -e "\nNumber of macro blocks:\t\t$nMacroBlocks\nnumber of rest blocks:\t\t$nRestBlocks"
	seqMacro=$(seq 1 $nMacroBlocks)
	seqRest=$(seq $(($lastBlock-$nRestBlocks+1)) $lastBlock)
	mkdir dock_result_blocksChunked_${sets[$setSel]}
	cd dock_result_blocksChunked_${sets[$setSel]}
	echo -e "Starting the merge\n"
	tot=0
	for i in $seqMacro
	do
		seqTmp=$(seq $((($i-1)*$numMacroC+1)) $(($i*$numMacroC)))
		for j in $seqTmp
		do
			cat ../${sets[$setSel]}/block$j/dock_result.sdf >> dock_result_macro$i.sdf
		done
		echo "macro block $i complete"
		tot=$(($tot+$(grep MOE dock_result_macro$i.sdf | wc -l | cut -d ' ' -f 1)))
		echo "tot counted: $tot"
	done
	for i in $seqRest
	do
		cat ../${sets[$setSel]}/block$i/dock_result.sdf >> dock_result_macro_rest.sdf
	done
	tot=$(($tot+$(grep MOE dock_result_macro_rest.sdf | wc -l | cut -d ' ' -f 1)))
	echo "tot DEFINITIVE counted: $tot"
	cd ..
fi


if [[ -d dock_result_blocksChunked_${sets[$setSel]} ]]
then
	merge=true
else
	merge=false
fi

cd $file_path


#activate env of DeepDock
source ~/envDeepDock/bin/activate
#EXTRACTION LABEL
python $DEEPDOCKNA/phase_2-3/Extract_labels.py -if False -n_it $currIt -protein $protein -file_path $file_path -t_pos $t_cpu -score S -zincid ZINCID -site $siteSel -set ${sets[$setSel]} -merge $merge


#if first iteration, then the dir iwill be the morgan dir

read -p "training, test and validation sets ready? [true|false] " ans
if $ans
then
	python $DEEPDOCKNA/phase_2-3/simple_job_models.py -n_it $currIt -mdd $morgan_directory -time 10:00:00 -file_path $file_path/$protein -nhp 24 -titr 11  -pfm 0.01 -plm 0.0001;
fi
