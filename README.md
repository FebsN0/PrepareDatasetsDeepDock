# PrepareDatasetsDeepDock

The following script are used to prepare training, testing and validation sets required by DeepDock-NonAutomated (forked and edited version available on my profile, https://github.com/FebsN0/Deep-Docking-NonAutomated.git)

NOTE!!!
The scripts suppose several things:
1) you have MOE (v2022.2 at least) installed on cluster with "Slurm Workload Manager" as job scheduler
2) you have prepared run.sh generated via "Batch" option inside Conformational Search MOE tool (example in csearchMainCode dir with a small example of input file)
3) you have prepared run.sh and rec.moe generated via "Batch" option inside Docking MOE tool (examples in dataXdockingStepBasedOnSite dir)
4) you have downloaded a database and properly prepared (washed, protonation and enantiomer state enumeration, etc). In my work I have downloaded the ZINC20 from https://files.docking.org/zinc20-ML/ (the database is already prepared!)
5) a log file with specific informations in a specific order (look in the example given further down)

CONSIDERATIONS:
Select option 1 or 2 for a generic cluster usage (options 3 and 4 are designed for specific cluster I have used which have internal problems, thus I have changed specific details for them)

EXAMPLE LOG FILE (delete the comment)
/home/fabiano/scratch/DeepDock                        #path main directory where smile, morganFP, projects, dataXdockingStepBasedOnSite and csearchMainCode are saved
gammaTubulin                                          #name of the project. It is actually just a directory, where everything will be saved
/home/fabiano/scratch/DeepDock/morganFP_original      #path of the directory where morganFP are saved (entire DB)
/home/fabiano/scratch/DeepDock/smileDir               #path of the directory where SMILE strings are saved (entire DB)
3000000                                               #number of size of all datasets (size of training + size of testing + size of validation sets)

Maybe in the future I will modify these scripts for a total generic usage, by using any conformational generator and docking program of one's choice
