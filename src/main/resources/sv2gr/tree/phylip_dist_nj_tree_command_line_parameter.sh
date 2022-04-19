#!/bin/bash

#must run as ./this.bash phylip_file working_dir output_dir distance_model outgroup_seq_index
if [ $# != 5 ]; then
	echo "./this	phylip_file working_dir output_dir distance_model[F84, Kimura, Jukes-Cantor, LogDet] outgroup_seq_index[int]"
	exit 1
fi

#full path to the phylip file that contains the multiple sequence alignment
#file name must end with .phy; 
#the base name of the file will be used to name all output files including the distance matric file and tree file
phylip_file=$1
file_basename=$(basename "$phylip_file" .phy) #file base name without path and extension 'phy'

#the directory where a temporary output directory will be created for the task
#should be pre-existing; if not , abort
working_dir=$2

#the output directory to store the resulted tree file for downstream analysis; should be different from the $working_dir
#if not exists, abort;
#this folder is normally used for storing of output files of multiple jobs that are run in the same batch
output_dir=$3


##################dnadist specific settings
distance_model=$4 #F84, Kimura, Jukes-Cantor, LogDet; F84 is the default and nothing needs to be added to command file; for Jukes-Contor, two 'D's should be put in the command file
#seq_is_sequential="true" #whether the phylip format contains sequential sequence or not (interleaving)
###################neighbor specific settings
outgroup_seq_index=$5 #index of the sequence used as outgroup root; if not given, the first one will be used as outgroup

###########
#the directory to change to when running phylip executables so that all generated command files and output files will be put there
#should be created by this bash and deleted after job is done
#note that phylip executables will automatically put the output files in the same directory from which they are run;
#thus, if there are multiple tasks to be run, each should be assigned a separate working directory to avoid overwriting on each other
tmp_out_dir_name=${file_basename}_tmp_out
tmp_out_dir=""

#the full path to the file contains all the commands when running the dnadist module interactively;
dnadist_command_file="" #phylip module 'dnadist'
neighbor_command_file="" #phylip module 'neighbor'

dnadist_matric_outfile="" #the automatically generated output file for calculated distance matrics
dnadist_matric_renamed_outfile=""

neighbor_outtree="" #the automatically generated nj tree file
neighbor_outfile="" #the automatically generated output file by neighbor
neighbor_renamed_outtree_file=""

#check if there is already running job in the working directory or not;
#note that due to the fact that phylip's output files are automatically generated with the same name, each job should be running in a separate working directory to avoid conflict!
preprocess(){
	#first check if phylip file exists; if not, abort
	if [ ! -f $phylip_file ]; then
		echo "Error: input phylip file $phylip_file does not exists!"
		exit 1
	fi
	
	
	#check if working and output directory exist 
	if [ ! -d $working_dir ]; then
		echo "Error: working directory $working_dir is not found!"
		exit 1
	fi
	
	if [ ! -d $output_dir ]; then
		echo "Error: output directory $output_dir is not found!"
		exit 1
	fi
	
	if [[ "$working_dir" == */ ]]; then #if working_dir ends with a '/'
		tmp_out_dir=${working_dir}${tmp_out_dir_name}/
	else
		tmp_out_dir=${working_dir}/${tmp_out_dir_name}/  #add a '/' after 
	fi
	
	#check if tmp out dir exists; if no, create it
	if [ ! -d $tmp_out_dir ]; then
		echo "creating temp out directory..."
		mkdir $tmp_out_dir
	fi
	
	dnadist_command_file=${tmp_out_dir}${file_basename}.dnadist_command_file #phylip module 'dnadist'
	neighbor_command_file=${tmp_out_dir}${file_basename}.neighbor_command_file #phylip module 'neighbor'

	dnadist_matric_outfile=${tmp_out_dir}outfile #the automatically generated output file for calculated distance matrics
	dnadist_matric_renamed_outfile=${tmp_out_dir}${file_basename}.${distance_model}.dist.matrics

	neighbor_outtree=${tmp_out_dir}outtree #the automatically generated nj tree file
	neighbor_outfile=${tmp_out_dir}outfile #the automatically generated output file by neighbor
	neighbor_renamed_outtree_file=${tmp_out_dir}${file_basename}.${distance_model}.dist.nj.tree.nwk

	#first check if command files already exists; if yes, abort
	if [ -f $dnadist_command_file ]; then
		echo "Error: dnadist command file already exists!"
		exit 1
	fi
	if [ -f $neighbor_command_file ]; then
		echo "Error: neighbor command file already exists!"
		exit 1
	fi
	
}

##create the command file for dnadist
create_dnadist_command_file(){
	
	#write to the dnadist command file
	##1. set up the input phylip file
	echo $phylip_file > $dnadist_command_file
	##2. set up distance_model
	if [ $distance_model = "Kimura" ]; then
		echo "D" >> $dnadist_command_file
	elif [ $distance_model = "Jukes-Cantor" ]; then
		echo $'D\nD' >> $dnadist_command_file #print two lines of 'D'
	elif [ $distance_model = "LogDet" ]; then
		echo $'D\nD\nD' >> $dnadist_command_file #print three lines of 'D'
	elif [ $distance_model != "F84" ]; then 
		echo "Error: unrecognized distance model $distance_model"
		exit 1
	fi
	##3. done
	echo "Y" >> $dnadist_command_file
}

run_phylip_dnadist(){
	ml PHYLIP/3.697-foss-2019b
	
	#phylip output files will be put in the same folder where the executable is run
	cd $tmp_out_dir
	
	#run dnadist
	dnadist < $dnadist_command_file
	
	#rename outfile
	mv $dnadist_matric_outfile $dnadist_matric_renamed_outfile
	
	#remove the command file
	rm $dnadist_command_file
	
	ml unload PHYLIP/3.697-foss-2019b
}

##create the command file for neighbor
create_neighbor_command_file(){
	#write to the dnadist command file
	##1. set up the input distance matrics file
	echo $dnadist_matric_renamed_outfile > $neighbor_command_file
	##2. set up outgroup_seq_index
	if [ $outgroup_seq_index != "" ]; then
		echo "O" >> $neighbor_command_file
		echo $outgroup_seq_index >> $neighbor_command_file
	fi
	
	##3. done
	echo "Y" >> $neighbor_command_file
}

run_phylip_neighbor(){
	ml PHYLIP/3.697-foss-2019b
	
	#phylip output files will be put in the same folder where the executable is run
	cd $tmp_out_dir
	
	#run dnadist
	neighbor < $neighbor_command_file
	
	#rename outfile
	mv $neighbor_outtree $neighbor_renamed_outtree_file
	
	#remove the command file
	rm $neighbor_command_file
	
	#remove the outfile
	rm $neighbor_outfile
	
	ml unload PHYLIP/3.697-foss-2019b
}


postprocess(){
	#move the distance file and tree file to the output dir
	mv $dnadist_matric_renamed_outfile $output_dir
	mv $neighbor_renamed_outtree_file $output_dir
	
	#remove the $tmp_out_dir
	cd $working_dir
	rm $tmp_out_dir -r
}


preprocess
create_dnadist_command_file
run_phylip_dnadist
create_neighbor_command_file
run_phylip_neighbor
postprocess