#!/usr/bin/env bash
#Set $genome_path
#Create genomes.list in genome_reference/ and write inside the names of the genomes line per line.
#Create "targets/target_n" as different combinations of pair samples you want to compare [sample1\tsample2], where "n" is the order of the genome in 'genomes.list'. If you want to compare the same combinations of samples for each genome, you must to create only one file "targets/target"
#$output_path	if you use various genomes, are going to create as many folders as reference genomes you use inside that folder.

. ~soft_bio_267/initializes/init_autoflow

output_path=$SCRATCH'/leng/new_det_sex'
actual_dir=`pwd`

 
while IFS= read reference; do

	output_workflow="`echo $output_path`/`echo $reference | cut -d "." -f 1`"
	
	echo -e "\nLaunching GetRegFlow for $reference\n"
 #Comment
	vars=`echo "
	\\$genome_file=$actual_dir/references/$reference,
	\\$samples_path=$actual_dir/samples,
	\\$samples=[male_1_aqg_pair;male_2_aqg_pair;female_1_aqg_pair;female_2_aqg_pair],
	\\$chunks=40,
	\\$path_target=$actual_dir/target,
	" | tr -d [:space:]`
	
	if [ ! -s $output_workflow ]; then
		mkdir -p $output_workflow
	fi
	
	if [ $1 == '1' ]; then	
	
		AutoFlow -w $actual_dir/mapping_comparison_template.af -o $output_workflow -m 16 -c 1 -s -V $vars $2
	
	elif [ $1 == '2' ]; then
	
		flow_logger -e $output_workflow -r all

	elif [ $1 == 'r' ]; then

		echo "Relaunching failed jobs"
		flow_logger -e $output_workflow -w -l

	fi

	#$actual_dir/mappingdd_comparison.sh 1 $genome $target $output_workflow $2
	#$actual_dir/mapping_comparison.sh 2 $genome $target $output_workflow
	#$actual_dir/mapping_comparison.sh r $genome $target $output_workflow

done < genomes.list
	
