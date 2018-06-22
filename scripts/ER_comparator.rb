#!/usr/bin/env ruby

require 'optparse'

CONTIG_NAME = 0
START = 1
STOP = 2 

#########################################################################################################
## METHODS
#########################################################################################################

def load_files(files, max_reads, program) #files is an array
	all_files = []
	if program == 'peakranger' 
		coverage_control_field = 8
		coverage_condition_field = 9
		files.each do |file|
			loaded_file = []
			File.open(file).each do |line|
				if line !~ /^#/
					line.chomp!
					if !line.empty?
						fields = line.split("\s")
						loaded_file << [fields[CONTIG_NAME], fields[START].to_i, fields[STOP].to_i] if fields[coverage_control_field].to_i <= max_reads || fields[coverage_condition_field].to_i <= max_reads
					end
				end
			end
			all_files << loaded_file
		end
	elsif program == 'ER_seeker'
		files.each do |file|
			loaded_file = []
			File.open(file).each do |line|
				line.chomp!
				if !line.empty?
					fields = line.split("\s")
					loaded_file << [fields[CONTIG_NAME], fields[START].to_i, fields[STOP].to_i]
				end
			end
			all_files << loaded_file
		end
	elsif program == 'epic'
		coverage_control_field = 3
		coverage_condition_field = 4
		files.each do |file|
			loaded_file = []
			line_number = 0
			File.open(file).each do |line|
				line.chomp!
				line_number += 1 
				if !line.empty? && line !~ /^#/ && line_number > 2
					fields = line.split("\s")
					loaded_file << [fields[CONTIG_NAME], fields[START].to_i, fields[STOP].to_i] if fields[coverage_control_field].to_i <= max_reads || fields[coverage_condition_field].to_i <= max_reads
				end
			end
			all_files << loaded_file
		end
	end
	return all_files
end


def find_coincidences(reference_record, current_record)
	record_to_save = nil
	if reference_record[CONTIG_NAME] == current_record[CONTIG_NAME]
		if current_record[START] >= reference_record[START] && 
			current_record[STOP] <= reference_record[STOP]
			record_to_save = current_record
		elsif current_record[START].between?(reference_record[START], reference_record[STOP]) 
			record_to_save = [reference_record[CONTIG_NAME], current_record[START], reference_record[STOP]]
		elsif current_record[STOP].between?(reference_record[START], reference_record[STOP])
			record_to_save = [reference_record[CONTIG_NAME], reference_record[START], current_record[STOP]]
		elsif current_record[START] < reference_record[START] && 
			current_record[STOP] > reference_record[STOP]
			record_to_save = reference_record
		end
	end
	return record_to_save
end


def compare_files(all_files) #all_files is a hash of arrays
	results_table = all_files.shift
	while !all_files.empty? do 
		current_table = all_files.shift
		records_to_save = []
		results_table.each do |saved_record|
			current_table.each do |current_record|
				record_to_save = find_coincidences(saved_record, current_record)
				records_to_save << record_to_save if !record_to_save.nil?
			end
		end
		results_table = records_to_save
	end
	return results_table
end


def save_results(filename, results_table)
	File.open("#{filename}", 'w') do |file|
		results_table.each do |record|
			file << "#{record.join("\t")}\n"
		end
	end
end




#########################################################################################################
## OPTS
#########################################################################################################
 
options = {}

OptionParser.new do  |opts|
	
	options[:input] = []
	opts.on("-i FILE", "--input", "Input file. Indicate between quotes") do |i|
		options[:input] = Dir.glob(i)
	end

	options[:program] = ''
	opts.on("-p COMPARATOR_PROGRAM_USED", "--program", "Set the used program for comparate BAM files. ARGUMENTS: 'epic' or 'peakranger'") do |p|
		options[:program] = p
	end

	options[:max_reads] = 0
	opts.on("-r INTEGER", "--reads", "Maximum reads that support one of the samples [DEFAULT = 0]") do |r|
		options[:max_reads] = r.to_i
	end
	options[:output] = nil
	opts.on("-o FILENAME", "--output", "Set output file") do |o|
		options[:output] = o
	end

	opts.on("-h", "--help", "Display options") do 
		puts opts
	end

end.parse!



#########################################################################################################
## MAIN
#########################################################################################################

abort("ERROR: The specified files not exist") if options[:input].empty?
abort("ERROR: The program used to comparate BAM files is not recognized") if options[:program] != 'epic' && options[:program] != 'peakranger' && options[:program] != 'ER_seeker'

all_files = load_files(options[:input], options[:max_reads], options[:program])
results_table = compare_files(all_files)
save_results(options[:output], results_table)
