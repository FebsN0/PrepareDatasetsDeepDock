#!/bin/bash

##############################################------------------------- FUNCTIONS -------------------------##############################################

function counterJobs(){
#$1 : increase rate
	countJobs=$(($countJobs+$1))
        if [[ $countJobs -ge $nRuns ]]
        then
        	echo -e "\n\nMAX PREPRARED JOBS based on nRUNS REACHED: exiting from the checking cycle\nLast Block checked:block$i\n"
                break
        fi
}

function checkHourMinSec (){
#check if the format is correct
	DD=$2
	HH=$(echo $1 | cut -d ':' -f 1)
	MM=$(echo $1 | cut -d ':' -f 2)
	SS=$(echo $1 | cut -d ':' -f 3)
	other=$(echo $1 | cut -d ':' -f 4)
	if [[ ! -n $HH || ! -n $MM || ! -n $SS || -n $other ]]
	then
		echo "ERROR: time wrong format"
             	exit 1
        elif [[ ! -n $DD && $HH -eq 0 && $MM -lt 20 ]]
	then
		echo "ERROR: too low minutes, at least 20"
                exit 1
	fi
}

function checkTime (){
#check day
	if [[ "$1" == *"-"* ]]
	then
		DD=$(echo $1 | cut -d '-' -f 1)
		if [[ $DD -gt 7 || $DD -lt 1 ]]
		then
			echo "ERROR: wrong format"
			exit 1
		fi
		HMS=$(echo $1 | cut -d '-' -f 2)
		checkHourMinSec $HMS $DD
	else
		checkHourMinSec $1
	fi
}

function selectSlurm (){
    if [[ -d confS/run.log ]]
    then
        slurmConfS=$(echo confS/run.log/*slurm/T1.err | awk '{print $NF}')
    elif [ -e confS/slurm* ]
    then
        slurmConfS=$(echo confS/slurm* | awk '{print $NF}')
    else
        slurmConfS=""
	fi
	if [[ -d dockG/run.log ]]
    then
        slurmDockG=$(echo dockG/run.log/*slurm/T1.err | awk '{print $NF}')
    elif [ -e dockG/slurm* ]
    then
        slurmDockG=$(echo dockG/slurm* | awk '{print $NF}')
    else
        slurmDockG=""
    fi
}

#find the name of block in which the job is still running. First line is the informations column
#specify if confSearch or Dock as input argument and type of job
# $1 : R or PD
# $2 : confSearch or Docking

function checkJobStatusBlock_i (){
	checkStatus=false
	echo -e "\tchecking the job status $2"
	sq | grep -P "^(?=.*$1)(?=.*$2)" > dataJobs.tmp
	case $1 in
		R)  awk '{print $1}' dataJobs.tmp > jobs.tmp;;
		PD) if [[ $cluster -eq 3 && $2 == "confSearch" ]]; then awk '{print $1}' dataJobs.tmp >  jobs.tmp; else awk '{print $1}' dataJobs.tmp | grep "\[1\]" | cut -d '_' -f 1 > jobs.tmp; fi;;
	esac
	for jobid in $(cat jobs.tmp)
	do
#sacct return different string: Narval as Beluga | Graham as Cedar
		case $cluster in
			1|3) checkinfo=$(sacct -j $jobid --format=workdir%100 | awk -F "/" '{printf "%s\t%s",$10,$9}');;
			2|4) checkinfo=$(sacct -j $jobid --format=workdir%100 | awk -F "/" '{printf "%s\t%s",$11,$10}');;
		esac
		if [[ $(echo $checkinfo | cut -d ' ' -f 1) == "block$i" && $(echo $checkinfo | cut -d ' ' -f 2) == "${sets[$setSel]}" ]]
                then
	 		echo -e "\tResult job status: $1"
			checkStatus=true
			break
                fi
	done
	if ! $checkStatus
	then
		echo -e "\tResult job status: not found"
	fi
	rm dataJobs.tmp jobs.tmp
}

# name array
# arrayBlocksReady2postprocess

function postprocess (){
#inside block$i
# $1 : from splitting part or normal part? [NORM for normal processing, otherwise none]
#CONVERSION MDB 2 SDF
#INPUT:         dock.mdb
#OUTPUT:        dock_result.sdf

#if normal processing or chuk part
	if [[ ! -n $1 ]]
    	then
        	cd dockG
        	moebatch -exec "db_ExportSD['dock.mdb','dock_result.sdf',['mol','ZINCID','S']]"
        	cd ..
        	mv dockG/dock_result.sdf .
    	else
		for ch in $allChunks
		do
			cd $ch
			moebatch -exec "db_ExportSD['dock.mdb','dock_result.sdf',['mol','ZINCID','S']]"
			cat dock_result.sdf >> $base/dock_result.sdf
			cd ..
        	done
		cd $base
    	fi
#find the missings ZINCID by comparing the initial file.smi (unluckly during the docking some compounds will be rejected) with the new file
    	awk '$2=="<ZINCID>"{ fetch_next = 1; next} fetch_next {print $1; fetch_next=0}' dock_result.sdf > block_res.zincid
    	mis=0; old=0; res=0;
#suppress second and third column, showing the unique ZINCID present in the first set. block_old.zincid was previously generated
 	comm -23 <(sort block_old.zincid) <(sort block_res.zincid) > missing.zincid
    	mis=$(wc -l missing.zincid | cut -d ' ' -f 1)
    	old=$(wc -l block_old.zincid | cut -d ' ' -f 1)
    	res=$(wc -l block_res.zincid | cut -d ' ' -f 1)
    	echo -e "\n\t\tCONCLUSION:\nmissing:\t$mis\ntotal:\t\t$old\nCompleted:\t$res" > result.log
#check if the results are okay before compressing
    	if [[ $res -eq 0 ]]
    	then
    		echo "$b : block_res.zincid = 0 ! ANOMALY" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
    	elif [[ $res -gt $old ]]
    	then
    		echo "$b : block_res.zincid > block_old.zincid MAGGIORE ANOMALY" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
    	elif [[ $res -le $(($old/100*98)) ]]
    	then
        	echo "$b : block_res.zincid < 98/100 block_old.zincid MINORE ANOMALY" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
    	else
#COMPRESS EVERYTHING if everything is fine. Sometime even csearch.mdb is ok, docking is done only to a really small amount of compounds. Something wrong!
		case $cluster in
			1|2) tar -zcf data.tar.gz missing* block* confS dockG --remove-files ;;
			3|4) tar -cf data.tar missing* block* confS dockG --remove-files ;;
        	esac
        	if [ -n $1  ]
        	then
            		echo -e "\t\t$b : postProcess SPLITTING completed" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
            		tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
        	else
            		echo "$b : postProcess completed"
        	fi
	fi
	count=$(($count-1))
}

# name array
# arrayBlocksReady2Docking
function docking (){
	cd dockG
	if [[ -d dockG/run.log ]]
	then
		run.log
	fi
	if [[ -e dock.mdb ]]
	then
		dock.mdb
	fi
	case $cluster in
        	1|2|4) sh run.sh -qsys slurm -qargs "--ntasks=10 --cpus-per-task=2 --mem=3G --nodes=1 --account=def-jtus --time=$timeDock --job-name=Docking" -submit;;
                3) sh run.sh -qsys slurm -qargs "--ntasks=1 --cpus-per-task=10 --mem=3G --nodes=1 --account=def-jtus --time=$timeDock --job-name=Docking" -submit;;
        esac
	count=$(($count-1))
}

# name array
# arrayBlocksReady2DockingResume
function resumeDocking(){
	cd dockG
        case $cluster in
                1|2|4) sh run.sh -qsys slurm -qargs "--ntasks=10 --cpus-per-task=2 --mem=3G --nodes=1 --account=def-jtus --time=$timeDockResume --job-name=DockingResume" -submit -- -resume;;
                3) sh run.sh -qsys slurm -qargs "--ntasks=1 --cpus-per-task=10 --mem=3G --nodes=1 --account=def-jtus --time=$timeDockResume --job-name=DockingResume" -submit -- -resume;;
        esac
	count=$(($count-1))
}

# name array
# arrayBlocksReadyConfSearchSPLITTED
#sometime, confSearch struggle to continue. Some the confSearch get stuck in calculation, thus, better split to save time

function restartConfSearchSPLIT(){
	currBlock=`pwd`
	cd confS
#remove old data
	case $cluster in
		1|2|4) rm -r run.log csearch.mdb;;
        3) rm -r slurm* csearch.mdb;;
	esac
	rows=$(wc -l ../block.smi | cut -d ' ' -f 1)
	sizeChunk=$(($rows/${nChunks[$idx]}))
	rest=$(($rows%$nChunks))
	echo ""
	echo -e "\n\tSPLITTING OF $b\nN elements:\t\t\t$rows\nN chunks:\t\t\t${nChunks[$idx]}\nSize of each smaller chunk:\t$sizeChunk\nrest:\t\t\t\t$rest" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
	tail -n5 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
#create subdir
	mkdir chunks
	cd chunks
#in the case there are residues
	if [[ $rest -ne 0 ]]
        then
		mkdir chunk_res
		tail -n $(($rows%${nChunks[$idx]})) ../../block.smi | tail -n $sizeChunk > chunk_res/chunk.smi
	fi
#prepare the splitted data
	sChunks=$(seq 1 ${nChunks[$idx]})
	for n in $sChunks
        do
                mkdir chunk_$n
                head -n $(($sizeChunk*$n)) ../../block.smi | tail -n $sizeChunk > chunk_$n/chunk.smi
	done
#start the slurm
	for ch in $(ls)
	do
		cd $ch
		moebatch -exec "db_Open['csearch_input.mdb','create']"
	    moebatch -exec "db_ImportASCII [ascii_file:'chunk.smi',db_file:'csearch_input.mdb',names:['mol','ZINCID'],types:['molecule','char'],titles:0]"
		case $cluster in
            1|2|4) cp ../../run.sh . ; sh run.sh -qsys slurm -qargs "--ntasks=10 --cpus-per-task=2 --mem=3G --nodes=1 --account=def-jtus --time=$timeConfS_SPLIT --job-name=confSearchRestartedSPLIT" -submit -- -i csearch_input.mdb;;
            3) cp ../../exeSbatchConfS.sh ../../run.sh . ; sed -i "s/#SBATCH --time=24:00:00/#SBATCH --time=$timeConfS_SPLIT/g" exeSbatchConfS.sh; sbatch exeSbatchConfS.sh;;
        esac
		cd ..
	done
	cd $currBlock
}


# name array
# arrayBlocksReady2ConfSearchRestart

function restartConfSearch(){
	cd confS
    case $cluster in
        1|2|4) rm -r run.log csearch.mdb; sh run.sh -qsys slurm -qargs "--ntasks=10 --cpus-per-task=2 --mem=3G --nodes=1 --account=def-jtus --time=$timeConfS --job-name=confSearchRestarted" -submit -- -i csearch_input.mdb;;
        3) rm -r slurm* csearch.mdb; sed -i 's/"#SBATCH --time=24:00:00"/"#SBATCH --time=$timeConfS"/g' exeSbatchConfS.sh; sbatch exeSbatchConfS.sh;;
	esac
	count=$(($count-1))
}

#function asking if restart the confSearch or splitting or none
function askRestartOrSplitConfSearch(){
    if [ -e confS/csearch_input.mdb ]
    then
        read -p "select: (1) restart the confSearch | (2) restart the confSearch splitted in smaller chunks | (other) none: " ans
        case $ans in
            1) arrayBlocksReady2ConfSearchRestart+=(block$i); echo -e "\tblock$i : confSearch restarted\n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log; counterJobs 1;;
            2) echo -e "\nblock$i has $(wc -l block.smi |cut -d ' ' -f 1) elements"; read -p "How many smaller chunks? " nChunks[$indexBlockSplit]; indexBlockSplit=$((indexBlockSplit+1)); arrayBlocksReadyConfSearchSPLITTED+=(block$i); echo -e "\tblock$i : confSearch SPLIT started\n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log; counterJobs ${nChunks[$(($indexBlockSplit-1))]};;
            *) echo -e "\tblock$i : confSearch not restarted\n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log;;
        esac
     else
        echo -e "\tblock$i : WTF.. csearch_input.mdb doesnt exist" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
        tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
    fi
}


function askRestartOrResumeDocking(){
    read -p "Restart the docking from start (1) or resume (2) or none (other)? " ans
    case $ans in
        1) arrayBlocksReady2Docking+=(block$i); echo -e "\tblock$i : Docking restarted\n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log; counterJobs 1;;
        2) arrayBlocksReady2DockingResume+=(block$i); echo -e "\tblock$i : Docking resumed \n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log; counterJobs 1;;
		*) echo -e "\tblock$i : Docking NO ACTION TAKEN \n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log;;
    esac
}




##############################################------------------------- INITIALIZATION -------------------------##############################################

#total is the total zincid including te not completed. totalProcessed is the total of only complete blocks
countJobs=0
totalold=0
missing=0
totalres=0
completed=0
nBlockProcessed=0
indexBlockSplit=0
idx=0

arrayBlocksReady2postprocess=()
arrayBlocksReady2Docking=()
arrayBlocksReady2DockingResume=()
arrayBlocksReady2ConfSearchRestart=()
arrayBlocksReadyConfSearchSPLITTED=()
nChunks=() #save the number of chunks for each block, its index is "indexBlockSplit"

##############################################------------------------- SETTINGS -------------------------##############################################

#select the cluster. In Graham and Narval all is fine. Cedar has problems in using more tasks (confSearch need moebatch function, docking need 1 task max).
#Moreover tar -cf/-xf on Cedar and Beluga instead of -zcf/-zxf on Graham and Narval (unknown reasons, gz compression not working)
echo -e "\n\t\t\tSETTINGs:\n"
echo -e "AVAILABLE CLUSTERs:\n 1 : Graham\n 2 : Narval\n 3 : Cedar\n 4 : Beluga"
read -p "select the cluster used currently: " cluster
case $cluster in
        1|2|3|4);;
        *) echo "wrong number. Select only 1, 2, 3 or 4"; exit 1;;
esac

# name project: ie name of directory where save anything
read -p 'nameDirectoryOfProject: ' protein
if [ ! -d $protein ]
then
        echo "the $protein directory doesn't exist! "
        exit 1
fi
#current iteration
read -p "current iteration: " currIt
if [ ! -d $protein/iteration_$currIt ]
then
        echo "the selected iteration directory doesnt exist"
        exit 1
fi
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

main=`pwd`
rm resultsStatus_site"$siteSel"_${sets[$setSel]}.log
cd $protein/iteration_$currIt/docking/site_$siteSel/${sets[$setSel]}

#parallelize some jobs, especially during the restart confSearch, starting/resuming docking jobs and postprocessing jobs
echo ""
read -p 'how many cpus use to parallelize? ' t_cpu

echo -e "\n"
#If at least 1 block exists, find the last processed block.
if [[ -n $(ls) ]]
then
        lastBlockProc=$(ls | sort -n -t 'k' -k 2 | tail -n 1)
        echo -e "the last block processed is:\t$lastBlockProc"
        lastRun="${lastBlockProc##*k}"          #the format is blockNUMBER
else
        lastBlockProc=0
        lastRun=0
        echo -e "\nThere are no processed blocks\nExiting because not possible to check and process. Run first the code 02\n"
        exit 0
fi

read -p "how many blocks runs? [max 1000 jobs] " nRuns
if [[ $nRuns -gt 1000 || $nRuns -lt 1 ]]
then
        exit 1
fi

#select the range of blocks to check
read -p "from what block start the process: " start
read -p "to what block end the process: " end

if [[ $end -gt $lastRun ]]
then
        end=$lastRun
fi
seqN=$(seq $start $end)


#check the number of elements of first and last block to figure out the time to sbatch (1000 ligs takes ~ 1 day)
if [[ -e block$start/data.tar.gz || -e block$start/data.tar ]]
then
        cd block$start
        case $cluster in
                1|2) tar -zxf data.tar.gz block.smi;;
                3|4) if [[ -e data.tar.gz ]]; then gunzip data.tar.gz; fi; tar -xf data.tar block.smi;;
        esac
        cd ..
fi

if [[ -e block$end/data.tar.gz || -e block$end/data.tar ]]
then
        cd block$end
        case $cluster in
                1|2) tar -zxf data.tar.gz block.smi;;
                3|4) if [[ -e data.tar.gz ]]; then gunzip data.tar.gz; fi; tar -xf data.tar block.smi;;
        esac
        cd ..
fi

startBlockNumElements=$(wc -l block$start/block.smi | cut -d ' ' -f 1)
endBlockNumElements=$(wc -l block$end/block.smi | cut -d ' ' -f 1)

if [[ -e block$start/data.tar.gz || -e block$start/data.tar ]]
then
        rm block$start/block.smi
fi
if [[ -e block$end/data.tar.gz || -e block$end/data.tar ]]
then
        rm block$end/block.smi
fi


echo -e "\nstart block has\t$startBlockNumElements\nend block has\t$endBlockNumElements\n\nChoose carefully the time"
read -p "time for ConfS (restart) [format: D-HH:MM:SS] " timeConfS
checkTime $timeConfS
read -p "time for ConfS (restart SPLITTED chunks) [format: D-HH:MM:SS] " timeConfS_SPLIT
checkTime $timeConfS_SPLIT
read -p "time for docking (start/restart) [format: D-HH:MM:SS] " timeDock
checkTime $timeDock
read -p "time for docking (resume) [format: D-HH:MM:SS] " timeDockResume
checkTime $timeDockResume
read -p "time for docking (SPLITTING) [format: D-HH:MM:SS] " timeDock_SPLIT
checkTime $timeDock_SPLIT
echo -e "\n\t\tSCRIPT STARTED\n\n"

##############################################------------------------- MAIN CODE -------------------------##############################################

pathBlocks=`pwd`
for i in $seqN
do
	nBlockProcessed=$(($nBlockProcessed+1))
#check if the block exists
	if [ ! -e block$i ]
    	then
        	echo -e "\tblock$i: ERROR1 : NOT PRESENT" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
        	tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
#if the block is not compressed (ie not complete), check what's the reason
    	elif [[ ! -e block$i/data.tar.gz && ! -e block$i/data.tar ]]
    	then
        	echo -e "block$i checking started"
        	cd block$i
        	base=`pwd`
        	totalold=$(($totalold+$(wc -l block.smi | cut -d ' ' -f 1)))
#exist directory of new slurm
		selectSlurm
#check if the block$i is not executed
        	if [[ ! -n $slurmConfS && ! -n $slurmDockG && ! -e confS/csearch.mdb && ! -e dockG/csearch.mdb && ! -e dockG/dock.mdb ]]
        	then
#check if chunks exist (if previously splitted). Check if completed, Cancelled, running, pending, complete

######################## CHUNK PART #######################à
#check if COnfSearch is done
            		if [[ -d confS/chunks && ! -d dockG/chunks ]]
            		then
                		complete=true
                		cd confS/chunks
                		allChunks=$(ls)
                		for ch in $allChunks
                		do
#if some of the chunk is not complete for any possible reason.
                    			case $cluster in
                        			1|2|4) slurmConfS_chunck=$(echo $ch/run.log/*slurm/T1.err | awk '{print $NF}');;
                        			3) slurmConfS_chunck=$(echo $ch/slurm* | awk '{print $NF}');;
                    			esac
                    			if [[ ! -n $(grep "ConfSearch *.* done" $slurmConfS_chunck) ]]
                    			then
                        			complete=false
                        			echo -e "\n##########################################################################\n"
                        			tail -n15 $slurmConfS_chunck
                        			echo -e "\n##########################################################################\n"
                        			echo -e "\t\tblock$i - $ch : ERROR16 : confSearch SPLITTING not complete" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                        			tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    			else
						echo ""
						echo -e "\t\tblock$i - $ch : confSearch SPLITTING complete" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                        			tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
					fi
				done
				cd $base
#if all chunks are complete, start the docking
				if $complete
				then
					mkdir dockG/chunks
					for ch in $allChunks
					do
						mkdir dockG/chunks/$ch
#copy run.sh and rec.moe. Move csearch.mdb. Start the job docking
						cp dockG/r* dockG/chunks/$ch
						mv confS/chunks/$ch/csearch.mdb dockG/chunks/$ch/.
						cd dockG/chunks/$ch
						case $cluster in
                            				1|2|4) sh run.sh -qsys slurm -qargs "--ntasks=10 --cpus-per-task=2 --mem=3G --nodes=1 --account=def-jtus --time=$timeDock_SPLIT --job-name=Docking_SPLIT" -submit;;
                					3) sh run.sh -qsys slurm -qargs "--ntasks=1 --cpus-per-task=10 --mem=3G --nodes=1 --account=def-jtus --time=$timeDock_SPLIT --job-name=Docking_SPLIT" -submit;;
                        			esac
						cd $base
					done
					echo ""
					echo -e "\t\tblock$i - $ch : Docking SPLITTING started" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    			tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
					counterJobs 4
				fi
#check if docking of chunks are complete
            		elif [[ -d confS/chunks && -d dockG/chunks ]]
            		then
                		complete=true
                		base=`pwd`
	        	        cd dockG/chunks
        	        	allChunks=$(ls)
                		for ch in $allChunks
                		do
#if some of the chunk is not complete for any possible reason.
                    			if [[ ! -n $(grep "Docking finished" $ch/run.log/*slurm/T1.err) ]]
                    			then
                        			complete=false
                        			echo -e "\n##########################################################################\n"
                        			tail -n50 $ch/run.log/*slurm/T1.err
                        			echo -e "\n##########################################################################\n"
                        			echo -e "\t\tblock$i - $ch : ERROR17 : Docking SPLITTING not complete" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                        			tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    			else
                        			echo ""
                        			echo -e "\t\tblock$i - $ch : Docking SPLITTING complete. PostProcess ready" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                        			tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    			fi
                		done
#postprocessing
                		if $complete
                		then
                    			b=block$i
                    			postprocess SPLIT
				fi
##################### END CHUNCK PART
#in cedar, slurm will be generated only after starting the job. Thus, "slurmConfS" is empty
            		elif [[ $cluster -eq 3 ]]
			then
				checkJobStatusBlock_i PD confSearch
				if $checkStatus
                		then
                    			echo "block$i : confSearch PENDING" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                		else
					echo -e "\tblock$i : ERROR2 : conformational search and docking not done" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    			tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    			askRestartOrSplitConfSearch
                		fi
			else
#GENERAL ERROR OF NOT STARTED
				echo -e "\tblock$i : ERROR2 : conformational search and docking not done" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
				tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
				askRestartOrSplitConfSearch
			fi

		######################### CONF_SEARCH PART ######################
#confSearch job exist, but not the docking job
		elif [[ -n $slurmConfS && ! -n $slurmDockG ]]
		then
#check if the block is under pending job or not. If not, errors occurred, thus required to be restarted
			if [[ ! -e confS/csearch.mdb && ! -e dockG/csearch.mdb ]]
			then
                if [[ -n $(grep "csearch_input.mdb is not readable" $slurmConfS) ]]
                then
                    echo -e "\tblock$i : ERROR3 : csearch_input un-readable. Required to re-create the input file. Restarted with normal settings" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    rm confS/csearch_input.mdb
                    moebatch -exec "db_Open['csearch_input.mdb','create']"
                    moebatch -exec "db_ImportASCII [ascii_file:'block.smi',db_file:'csearch_input.mdb',names:['mol','ZINCID'],types:['molecule','char'],titles:0]"
                    mv csearch_input.mdb confS/.
                    arrayBlocksReady2ConfSearchRestart+=(block$i)
                    counterJobs 1
                else
                    checkJobStatusBlock_i PD confSearch
                    if $checkStatus
                    then
                        echo "block$i : confSearch PENDING" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    else
                        tail -n50 $slurmConfS
						echo -e "\n##########################################################################\n"
                        echo -e "\tblock$i : ERROR4 : csearch.mdb doesnt exist, the job is terminated for unknown reasons" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                        tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                        askRestartOrSplitConfSearch
                    fi
				fi
#check if the job was cancelled. Maybe x time expired and other errors. Other errors could occurs. In case of choosing SPLITTED
            elif [[ -e confS/csearch.mdb && -n $(grep "CANCELLED" $slurmConfS) ]]
			then
                tail -n50 $slurmConfS
				echo -e "\n##########################################################################\n"
                echo -e "\tblock$i : ERROR5 : csearch job CANCELLED" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
				tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                askRestartOrSplitConfSearch

#check if confSearch is running. Sometime happens that the job is terminated without reporting any error. No absolute reason reported. Just finished.
#the slurm log is like that the job is still running, but actually it is not true (no any job in squeue)
			elif [[ -e confS/csearch.mdb && ! -n $(grep "ConfSearch *.* done" $slurmConfS) ]]
			then
                checkJobStatusBlock_i R confSearch
				if $checkStatus
                then
                    echo "block$i : confSearch is running" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                else
                    tail -n20 $slurmConfS
					echo -e "\n##########################################################################\n"
                    echo -e "\tblock$i : ERROR6: csearch.mdb exist, the job is terminated for some reasons" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
					askRestartOrSplitConfSearch
				fi
#confSearch is complete
			elif [[ -e confS/csearch.mdb && -n $(grep "ConfSearch *.* done" $slurmConfS) ]]
            then
#check if the number of molecule is greater than 98% by observing slurm job
                numLigsConfS_processed=$(grep "Processed" $slurmConfS | tail -n1 | awk '{print $2}')
                old=$(wc -l block_old.zincid | cut -d ' ' -f 1)
                if [[ $numLigsConfS_processed -le $(($old/100*98)) ]]
                then
                    tail -n20 $slurmConfS
                    echo -e "\n##########################################################################\n"
                    echo -e "\tblock$i : ERROR7 confS returned num ligands < 98/100 of block_old.zincid ANOMALY" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    askRestartOrSplitConfSearch
                else
#all is good
                    mv confS/csearch.mdb dockG/.
                    arrayBlocksReady2Docking+=(block$i)
                    counterJobs 1
                fi
                    #csearch.mdb already moved, but for some reasons, the docking job dir is not created. Maybe PENDING. Infact, during pending, run.log is created, but not T1.err, thus $slurmDockG is empty
			elif [[ -e dockG/csearch.mdb && ! -e dockG/dock.mdb ]]
			then
				checkJobStatusBlock_i PD Docking
				if $checkStatus #in pending status, slurm inside run.log doesnt exist yet
                then
                    echo "block$i : Docking PENDING" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                else
                    echo -e "\tblock$i : ERROR8 : csearch.mdb is moved already, but the docking job has not started for some reasons. Docking restarted" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
					tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    arrayBlocksReady2Docking+=(block$i)
                    counterJobs 1
				fi
			else
#some unpredicted error occurs
				tail -n50 $slurmConfS
                echo -e "\n##########################################################################\n"
				echo -e "\tblock$i : ERROR9 unknown" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
				tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
				echo -e "\t\tInside block$i:"
				ls
				echo -e "\t\tInside confS dir:"
				ls confS
				echo -e "\t\tInside dockG dir:"
				ls dockG
                read -p "ERROR9 - Restart the confSearch (1) or SPLIT the confSearch (2) or the docking (3) or none (other)? " ans
                case $ans in
                    	1) arrayBlocksReady2ConfSearchRestart+=(block$i); echo -e "\tblock$i : confSearch restarted\n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log; counterJobs 1;;
                    	2) echo -e "\nblock$i has $(wc -l block.smi |cut -d ' ' -f 1) elements"; read -p "How many smaller chunks? " nChunks[$indexBlockSplit]; indexBlockSplit=$((indexBlockSplit+1)); arrayBlocksReadyConfSearchSPLITTED+=(block$i); echo -e "\tblock$i : confSearch SPLIT started\n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log; counterJobs ${nChunks[$(($indexBlockSplit-1))]};;
			3) arrayBlocksReady2Docking+=(block$i); echo -e "\tblock$i : Docking restarted\n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log; counterJobs 1;;
			*) echo -e "\t\tblock$i : NO ACTION TAKEN \n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log;;
				esac
			fi
		######################### DOCKING PART ######################
#if the job of docking is somehow done
		elif [[ -n $slurmDockG ]]
        then
            if [[ -e dockG/csearch.mdb && ! -e dockG/dock.mdb ]]
			then
                checkJobStatusBlock_i PD Docking
#check if the block is under pending docking job or not. If first time, the result is always false. It is true in the previous call (confoSearch section)
				if $checkStatus
                then
                    echo "block$i : Docking PENDING" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
#in some cluster, MPU dont start for some reasom, thus the job fails
                elif [[ -n $(grep "MPU initialization failed" $slurmDockG) || -n $(grep "MPU multi-processor was shut down" $slurmDockG) ]]
				then
                    echo -e "\tblock$i : ERROR10: docking job failed to start because of MPU initialization fail or MPU multi-processor was shut down. Restarted" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
					arrayBlocksReady2Docking+=(block$i)
					counterJobs 1
				else
					echo -e "\n##########################################################################\n"
					tail -n50 $slurmDockG
                    echo -e "\n##########################################################################\n" 
					echo -e "\tblock$i : ERROR11: csearch.mdb is in dockG dir, but the docking job has not started for some reasons." >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
					askRestartOrResumeDocking
                fi
#dock.mdb exist, but check if docking job is cancelled, maybe because time expired
			elif [[ -e dockG/dock.mdb && -n $(grep "CANCELLED" $slurmDockG) ]]
			then
				echo -e "\n##########################################################################\n"
				grep Ligands $slurmDockG | tail
                echo -e "\n##########################################################################\n"
				echo -e "\tblock$i : ERROR12: docking job CANCELLED" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
				askRestartOrResumeDocking
#if docking is running
			elif [[ -e dockG/dock.mdb && ! -n $(grep "Docking finished" $slurmDockG) && -n $(grep "Docking started" $slurmDockG) ]]
			then
				checkJobStatusBlock_i R Docking
				if $checkStatus
                then
                    echo "block$i : Docking still running" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
				else
#sometime jobs stops without reporting any reasons, the slurm is like it's continuing to process but the job is not actually running
                    echo -e "\n##########################################################################\n"
                    tail -n50 $slurmDockG
                    echo -e "\n##########################################################################\n"
					echo -e "\tblock$i : ERROR13 : dock.mdb exists, but the docking job has stopped for unknown reasons and CANCELLED is not reported." >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
                    tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
					askRestartOrResumeDocking
                fi
#if docking is complete. ready to postprocess
            elif [[ -n $(grep "Docking finished" $slurmDockG) && -e dockG/dock.mdb && ! -e result.log ]]
            then
                echo "block$i : docking complete, ready to post-process"
		arrayBlocksReady2postprocess+=(block$i)
#sometimes dock.mdb is created even if inside run.log there is no Docking started
            else
                echo -e "\tblock$i : ERROR14 the job gets stuck for unknown reasons. " >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
				tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
				askRestartOrResumeDocking
		fi
		else
			echo -e "\tblock$i : ERROR15 unknown : problems with $slurmConfS and $slurmDockG" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
			tail -n1 $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
		fi

#END CHECK
		cd $pathBlocks
		echo -e "block$i checking finished\n"

#EVERYTHING SHOULD BE DONE and FINE. ConfSearch and Docking are complete (thereotically). compressed files exist
	elif [[ -e block$i/data.tar.gz || -e block$i/data.tar ]]
	then
        	cd block$i
        #START CANCEL PART
#PARTE DA CANCELLARE DOPO. POICHE CHECKING VIENE FATTO NELLA FNZIONE postprocessing
		if [ ! -e result.log ]
		then
			case $cluster in
				1|2) tar -zxf data.tar.gz missing.zincid block.smi block_res.zincid ;;
				3|4) if [[ -e data.tar.gz ]]; then gunzip data.tar.gz; fi; tar -xf data.tar missing.zincid block.smi block_res.zincid;;
			esac
#extract info
			old=$(wc -l block.smi | cut -d ' ' -f 1)
			res=$(wc -l block_res.zincid | cut -d ' ' -f 1)
			mis=$(wc -l missing.zincid | cut -d ' ' -f 1)
			missing=$(($missing+$mis))
			echo -e "\n\t\tCONCLUSION:\nmissing:\t$mis\ntotal:\t\t$old\nCompleted:\t$res" > result.log
			rm block* missing*
		else
		#END CANCEL PART
			old=$(grep total result.log | awk '{print $2}')
			res=$(grep Completed result.log | awk '{print $2}')
			missing=$(($missing + $(grep missing result.log | awk '{print $2}')))
		fi
		if [[ $res -eq 0 ]]
		then
			echo "block$i : block_res.zincid = 0 ! ANOMALY" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
		elif [[ $res -gt $old ]]
		then
			echo "block$i : block_res.zincid > block_old.zincid MAGGIORE ANOMALY" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
		elif [[ $res -le $(($old/100*98)) ]]
        	then
            	echo "block$i : block_res.zincid < 98/100 block_old.zincid MINORE ANOMALY" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
		else
#EVERYTHING WORKED WELL......
			totalold=$(($totalold+$old))
			totalres=$(($totalres+$res))
			completed=$(($completed+1))
			if [[ ! -e dock_result.sdf ]]
			then
				case $cluster in
                        	        1|2) tar -zxf data.tar.gz dockG/dock_result.sdf;;
                                	3|4) if [ -e data.tar.gz ]; then gunzip data.tar.gz; fi; tar -xf data.tar dockG/dock_result.sdf;;
				esac
				mv dockG/dock_result.sdf .
			fi
			echo "block$i OK"
		fi
		cd ..
	else
		echo "\nblock$i : ERROR18 : reason unknown"
	fi
done

########################################### START THE PROCESSING (confSearch-Docking-Postprocess based on built arrays)#############################

echo -e "\n\t\tSTART THE PROCESS (confSearch/Docking/postProcess)\n"
cd $pathBlocks

# POSTPROCESS
if [[ -n ${arrayBlocksReady2postprocess[@]} ]]
then
    echo -e "blocks X postprocessing:\n\t ${arrayBlocksReady2postprocess[@]}"
    read -p "start the postprocess? [y|other]: " ans
    if [[ $ans == y ]]
    then
        echo ""
        count=1
        for b in ${arrayBlocksReady2postprocess[@]}
        do
            cd $b
            postprocess &
            if [[ $(($count%$t_cpu)) -eq 0 ]]
            then
                wait
            fi
            count=$(($count+1))
            cd $pathBlocks
        done
        wait
    fi
fi

#CONF SEARCH - RESTART
if [[ -n ${arrayBlocksReady2ConfSearchRestart[@]} ]]
then
    echo -e "\nblocks X confSearch:\n\t${arrayBlocksReady2ConfSearchRestart[@]}"
    read -p "start the conf search? [y|other]: " ans
    if [[ $ans == y ]]
    then
        echo ""
        count=1
        for b in ${arrayBlocksReady2ConfSearchRestart[@]}
        do
            cd $b
            restartConfSearch &
            if [[ $(($count%$t_cpu)) -eq 0 ]]
            then
                wait
            fi
            count=$(($count+1))
            cd $pathBlocks
        done
        wait
    fi
fi

# CONF SEARCH - SPLITTING. No parallelization because there is manual preparation
if [[ -n ${arrayBlocksReadyConfSearchSPLITTED[@]} ]]
then
    echo -e "\nblocks X confSearch SPLIT:\n\t${arrayBlocksReadyConfSearchSPLITTED[@]}"
    read -p "start the conf search - SPLIT? [y|other]; " ans
    if [[ $ans == y ]]
    then
        echo "array nChunks: ${nChunks[@]}"
        for b in ${arrayBlocksReadyConfSearchSPLITTED[@]}
        do
            cd $b
            restartConfSearchSPLIT
            cd $pathBlocks
            idx=$(($idx+1))
        done
    fi
fi

# DOCKING - START/RESTART
if [[ -n ${arrayBlocksReady2Docking[@]} ]]
then
    echo -e "\nblocks X docking:\n\t${arrayBlocksReady2Docking[@]}"
    read -p "start the docking? [y|other]: " ans
    if [[ $ans == y ]]
    then
    count=1
    for b in ${arrayBlocksReady2Docking[@]}
    do
        cd $b
        docking &
        if [[ $(($count%$t_cpu)) -eq 0 ]]
        then
            wait
        fi
        count=$(($count+1))
        cd $pathBlocks
    done
    wait
    fi
fi

# DOCKING - RESUME
if [[ -n ${arrayBlocksReady2DockingResume[@]} ]]
then
    echo -e "\nblocks X docking - resume:\n\t${arrayBlocksReady2DockingResume[@]}"
    read -p "start the resume of docking? [y|other]: " ans
    if [[ $ans == y ]]
    then
        count=1
        for b in ${arrayBlocksReady2DockingResume[@]}
        do
            cd $b
            resumeDocking &
            if [[ $(($count%$t_cpu)) -eq 0 ]]
            then
                wait
            fi
            count=$(($count+1))
            cd $pathBlocks
        done
        wait
    fi
fi

echo -e "\n\t\tFINISHED\n"

echo -e "\n\tCONCLUSION CHECKING\n" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
echo -e "Total block_old.zincid:\t\t\t\t$totalold\nTotal block_res.zincid:\t\t\t\t$totalres\t\t(NB: include only blocks with no problems)\nTotal missing:\t\t\t\t\t$missing" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log
echo -e "Blocks completed over $nBlockProcessed blocks processed:\t\t\t$completed" >> $main/resultsStatus_site"$siteSel"_${sets[$setSel]}.log

