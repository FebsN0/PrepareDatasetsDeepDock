#!/bin/bash

#############-----------------------------------PREPARATION X DOCKING + LAUNCH SBATCH--------------------------#####################

#since ComputeCanada dont allow the generation of more than 1 million of files, needed to split in different docking blocks
#each block will be processesed indipendentely. Docking with MOE

###NOTE!!!!!!!!!!!!!!!!
# moebatch -mpu is not working properly. Since it is necessary to use sh run.sh
# this script start the slurm job of CONFORMATIONAL SEARCH ONLY! 


############---------------------- FUNCTIONs----------------------############################

function extract_convert_confSearch(){
#extract N chosen ligs from dataset generated in previous step phase1_A. Save it in block$j
#check if it's the last block, which represent residue. Not extract 10 anymore
#every time this program start, count start from 1.
	if [ $j -eq $(($lastRun+$nBlocks+1)) ]
	then
		tail -n $residue $smileDir/${sets[setSel]}_smiles_final_updated.smi > block.smi
	else
		head -n $(($totalProcessed+($nElemXblock*$count))) $smileDir/${sets[setSel]}_smiles_final_updated.smi | tail -n $nElemXblock > block.smi
	fi
	cut -d ' ' -f 2 block.smi > block_old.zincid
#conversion SMILE 2 MDB with db_ImportASCII
        moebatch -exec "db_Open['csearch_input.mdb','create']"
        moebatch -exec "db_ImportASCII [ascii_file:'block.smi',db_file:'csearch_input.mdb',names:['mol','ZINCID'],types:['molecule','char'],titles:0]"

#prepare data and copy the main code of conformational search into subdirectory
	mkdir confS
	cd confS
	mv ../csearch_input.mdb .
	cp $file_path/csearchMainCode/run.sh .
	cd ..

#prepare data X docking: copy from main dir rec.moe and run.sh based on selected site
#the input for the docking will be generated during SBATCH process
	mkdir dockG
	cd dockG
        cp $file_path/dataXdockingStepBasedOnSite/dockSite$siteSel/* .
#START WITH CONFORMATIONAL SEARCH
	cd ../confS
	case $cluster in
		3) cp $file_path/sample_exeConfSearch.sh exeSbatchConfS.sh; sed -i 's/"#SBATCH --time=24:00:00"/"#SBATCH --time=$timeConfS"/g' exeSbatchConfS.sh; sbatch exeSbatchConfS.sh;;
		1|2|4) sh run.sh -qsys slurm -qargs "--ntasks=10 --cpus-per-task=2 --mem-per-cpu=2G --nodes=1 --account=def-jtus --time=$timeConfS --job-name=confSearch" -submit -- -i csearch_input.mdb
	esac
	count=$(($count-1))
	echo "block$j : preparation completed"
}

#############----------------------------------SETTINGs --------------------------####smileDir#################

# $1: file logs.txt
if [[ ! -n $1 ]]
then
        echo "logs.txt missing"
        exit 1
fi

echo -e "\n\t\t\tSETTINGs:\n"
echo -e "AVAILABLE CLUSTERs:\n 1 : Graham\n 2 : Narval\n 3 : Cedar\n 4 : Beluga"
read -p "select the cluster used currently: " cluster
case $cluster in  
        1|2|3|4);;
        *) echo "wrong number. Select only 1, 2, 3 or 4"; exit 1;;
esac
echo ""

# filepath definition
sets=(train test valid)
file_path=`sed -n '1p' $1`

# name project: ie name of directory where save anything
protein=`sed -n '2p' $1`
read -p "name project is $protein, change or not? [0: change | other/none: ok] " ans
case $ans in
        0) read -p 'nameDirectoryOfProject: ' protein;;
        *) ;;
esac
if [ ! -d $file_path/$protein ]; then echo "$protein directory doesn't exist!"; fi


read -p "current iteration: " currIt
if [ ! -d $file_path/$protein/iteration_$currIt ]
then
        echo "Iteration $currIt dir inside $protein dir doesnt exist"
        exit 1
fi

#name smile dir
smileDir=$file_path/$protein/iteration_$currIt/smile

#select the binding site found previously
echo -e "\nBINDING SITES AVAILABLE:\n 1 (data saved on graham)\n 2 (data saved on narval)\n 3 (data saved on cedar)\n 4 (data saved on beluga)"
read -p 'choose the binding site? ' siteSel
case $siteSel in
        1|2|3|4);;
        *) echo "wrong number. Select only 1,2,3 or 4"; exit 1;;
esac

echo -e ""
sets=(train test valid)
echo -e "DATASETs:\n 0 (TRAINING SET - train Dir)\n 1 (TEST SET - test Dir)\n 2 (VALIDATION SET - valid Dir)"
read -p 'which dataset to run? ' setSel
case $setSel in
        0|1|2);;
        *) echo "wrong answer. Select only 0,1 or 2"; exit 1;;
esac

echo ""
read -p 'number of cores to parallelize this script: ' t_cpu
echo ""

cd $file_path/$protein/iteration_$currIt
if [ ! -d docking ]
then
	mkdir docking
	mkdir docking/site_$siteSel
	mkdir docking/site_$siteSel/${sets[$setSel]}
elif [ ! -d docking/site_$siteSel ]
then
	mkdir docking/site_$siteSel
	mkdir docking/site_$siteSel/${sets[$setSel]}
elif [ ! -d docking/site_$siteSel/${sets[$setSel]} ]
then
	mkdir docking/site_$siteSel/${sets[$setSel]}
fi
cd docking/site_$siteSel/${sets[$setSel]}


#if already saved
totalLigands=$(wc -l $smileDir/${sets[setSel]}_smiles_final_updated.smi | cut -d ' ' -f 1)
#count how many ligands are processed in the past. In the case of accidental elimination of orginal file
totalProcessed=0
if [[ -e ../totalProcessed_${sets[setSel]}.log && -n $(ls) ]]
then
#last row contain the last updated number of extracted ligands
	totalProcessed=$(tail -n1 ../totalProcessed_${sets[setSel]}.log)
#it was never counted  before. Create new file log
elif [[ -n $(ls) ]]
then
	for b in block*
	do
		cd $b
#if a process is still running during the counting, it means block_old is not compressed yet. Otherwise, it is in the data.tar.gz
		if [[ -e data.tar.gz || -e data.tar ]]
		then
			case $cluster in
				1|2) tar -zxf data.tar.gz block,smi;;
                                3|4) tar -xf data.tar block,smi;;
                        esac
			totalProcessed=$(($totalProcessed+$(wc -l block.smi | cut -d ' ' -f 1)))
			rm block.smi
		elif [[ -e block.smi ]]
		then
			totalProcessed=$(($totalProcessed+$(wc -l block.smi | cut -d ' ' -f 1)))
		else
			echo "something is wrong in $b ..."
			exit 1
		fi
		cd ..
	done
	echo -e "total ligands processed previously:\n$totalProcessed" > ../totalProcessed_${sets[setSel]}.log
fi

if [[ $totalLigands -le $totalProcessed ]]
then
	echo "all possible ligands in site $siteSel and in ${sets[$setSel]} set are already processed! "
        exit 0
else
	echo "Total number of processed ligands previously: $totalProcessed"
fi


read -p 'number of ligands for each blocks? ' nElemXblock
echo -e "\nChoose carefully the time based on number of ligands (1000 ligs/day)"
read -p "time for ConfS (restart) [format: D-HH:MM:SS] " timeConfS

#generate last block, count them
residue=$((($totalLigands-$totalProcessed)%$nElemXblock))
nBlocks=$((($totalLigands-$totalProcessed)/$nElemXblock))
#there are two types of blocks: mostly contains 10 elements. The last block depend on residue elements

echo -e "\n\ntotal in\t\t\t ${sets[$setSel]}_smiles_final_updated.smi:\t\t$totalLigands" >> ../splittingSmile_${sets[$setSel]}.log
echo -e "total extracted from\t\t ${sets[$setSel]}_smiles_final_updated.smi:\t\t$totalProcessed" >> ../splittingSmile_${sets[$setSel]}.log
echo -e "residues in\t\t\t ${sets[$setSel]}_smiles_final_updated.smi:\t\t$residue" >> ../splittingSmile_${sets[$setSel]}.log
echo -e "nBlocks with $nElemXblock elements in\t ${sets[$setSel]}_smiles_final_updated.smi:\t\t$nBlocks" >> ../splittingSmile_${sets[$setSel]}.log
cat ../splittingSmile_${sets[setSel]}.log
echo -e "\n"

#if at least 1 block exist
if [[ -n $(ls) ]]
then
#find the last processed block
        lastBlockProc=$(ls | sort -n -t 'k' -k 2 | tail -n 1)
        echo "the last block processed is: $lastBlockProc"
       	lastRun="${lastBlockProc##*k}"          #the format is blockNUMBER
else
	lastBlockProc=0
	lastRun=0
        echo -e "\nthere are no processed blocks! \nthus, lastRun=$lastRun\n"
fi

read -p "how many blocks runs? [max 1000 jobs, consider that each run processes $nElemXblock ligs ] " nRuns
if [[ $nRuns -gt 1000 || $nRuns -lt 1 ]]
then
        exit 1
fi

#prepare array of blocks to be created based on the last created block and the number of nRuns. Start must not overcome the last possible block+1
if [[ $(($nRuns)) -gt $(($nBlocks+1)) ]]
then
	end=$(($lastRun+$nBlocks+1))
else
	end=$(($lastRun+$nRuns))
fi

if [[ $(($lastRun+1)) -gt $end ]]
then
	echo -e "\n\tMAX NUMBER OF BLOCKS REACHED. EVERY LIGAND IS ALREADY PROCESSED"
	exit 0
else
	start=$(($lastRun+1))
fi

seqBlocks=$(seq $start $end)
echo -e "\n\tSTART BLOCK:\t$start\n\tEND BLOCK:\t$end\n"
read -p "setting completed, continue main code? [yes|no]" answer
if [[ $answer != yes ]]
then
	echo -e "REJECTED\n\n" >> ../splittingSmile_${sets[setSel]}.log
	exit 0
else
	echo -e "STARTED\n\n" >> ../splittingSmile_${sets[setSel]}.log
fi
echo -e "\nSTART\n"

###########################----------------------- MAIN CODE -----------------##################

blockDir=`pwd`
count=1
for j in $seqBlocks
do
	mkdir block$j
        cd block$j
#extract the 10 ligs, called as tmp.smi
	extract_convert_confSearch &
        if [ $(($count%$t_cpu)) -eq 0 ]
        then
        	wait
        fi
        count=$(($count+1))
	cd $blockDir
done
wait



#update total number processed ligands
for j in $seqBlocks
do
	totalProcessed=$(($totalProcessed+$(wc -l block$j/block.smi | cut -d ' ' -f 1)))
done
echo $totalProcessed >> ../totalProcessed_${sets[setSel]}.log
