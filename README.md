# PrepareDatasetsDeepDock

The following script are used to prepare training, testing and validation sets required by DeepDock-NonAutomated (forked and edited version available on my profile, https://github.com/FebsN0/Deep-Docking-NonAutomated.git)

NOTE!!!
The scripts suppose several things:
1) you have MOE (v2022.2 at least) installed on cluster with "Slurm Workload Manager" as job scheduler
2) you have prepared run.sh generated via "Batch" option inside Conformational Search MOE tool (example in csearchMainCode dir with a small example of input file)
3) you have prepared run.sh and rec.moe generated via "Batch" option inside Docking MOE tool (examples in dataXdockingStepBasedOnSite dir)
4) a log file (see an example on logs.txt)
5) you have downloaded a database and properly prepared (washed, protonation and enantiomer state enumeration, etc). In my work I have downloaded the ZINC20 from https://files.docking.org/zinc20-ML/ (the database is already prepared!)

Select option 1 or 2 for a generic cluster usage (options 3 and 4 are designed for specific cluster I have used which have internal problems, thus I have changed specific details for them)

Maybe in the future I will modify these scripts for a total generic usage, by using any conformational generator and docking program of one's choice
