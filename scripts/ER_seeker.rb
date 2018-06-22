#!/usr/bin/env ruby
#The program takes two sorted depth_per_nucleotide files from `samtools sort | samtools depth -a`, it divide each file in mapped and non-mapped regions, and it compare all the files for extract enriched regions with no mapping in one of them 
require 'optparse'
COVERAGE = 1
COORDINATE = 0
START = 0
STOP = 1
TAG = 2
######################################################################################################################
## METHODS
#############################################################################################################################
def load_file(filename, min_coverage) #files is an array with the filenames
	chrs = {}
	regions = []
	reg = []
	last_chr = nil
	last_coordinate = nil
	last_mapped = nil
	File.open(filename).each do |line|
		line.chomp!
		fields = line.split("\t")
		current_chr = fields.shift
		coordinate = fields[COORDINATE].to_i
		coverage = fields[COVERAGE].to_i
		if coverage > min_coverage
			mapped = true
		else
			mapped = false
		end
		if last_chr != current_chr && !last_chr.nil?
			regions << reg.concat([last_coordinate, last_mapped])
			chrs[last_chr] = regions 
			reg = []
			regions = []
		elsif mapped != last_mapped && !last_mapped.nil?
			regions << reg.concat([last_coordinate, last_mapped])
			reg = []
		end
		reg << coordinate if reg.empty?
		last_chr = current_chr
		last_coordinate = coordinate
		last_mapped = mapped
	end
	regions << reg.concat([last_coordinate, last_mapped])
	chrs[last_chr] = regions
	return chrs
end

def compare_samples(control, treatment, min_length)
	enriched_chrs = {}
	control.each do |control_chr, control_regions|
		enriched_regions = []
		treatment_regions = treatment[control_chr]	
		control_regions.each do |control_region|
			treatment_regions.each do |treatment_region|
				enriched_region = compare_records(control_region, treatment_region, min_length)
				enriched_regions << enriched_region if !enriched_region.nil?
			end
		end
		enriched_chrs[control_chr] = enriched_regions if !enriched_regions.empty?
	end
	return enriched_chrs
end

def compare_records(control_record, treatment_record, min_length)
	record_to_save = nil
	if control_record[TAG] != treatment_record[TAG]
		if treatment_record[START] >= control_record[START] && 
			treatment_record[STOP] <= control_record[STOP]
			record_to_save = [treatment_record[START], treatment_record[STOP]]
		elsif treatment_record[START].between?(control_record[START], control_record[STOP]) 
			record_to_save = [treatment_record[START], control_record[STOP]]
		elsif treatment_record[STOP].between?(control_record[START], control_record[STOP])
			record_to_save = [control_record[START], treatment_record[STOP]]
		elsif treatment_record[START] < control_record[START] && 
			treatment_record[STOP] > control_record[STOP]
			record_to_save = [control_record[START], control_record[STOP]]
		end
		if !record_to_save.nil? && record_to_save.last - record_to_save.first < min_length
			record_to_save = nil
		end
	end
	return record_to_save
end

def save_output(enriched_chrs)
	enriched_chrs.each do |chr_name, enriched_regions|
		enriched_regions.each do |enriched_region| 
			puts "#{chr_name}\t#{enriched_region.join("\t")}"
		end
	end
end
#############################################################################################################################
## INPUT PARSING
#############################################################################################################################
options = {}

OptionParser.new do  |opts|

	options[:control] = ''
	opts.on("-c control_FILE", "--file", "Set the file to extract info.") do |c|
		options[:control] = c
	end

	options[:treatment] = ''
	opts.on("-t treatment_FILE", "--file", "Set the file to extract info.") do |t|
		options[:treatment] = t	
	end 

	options[:min_coverage] = 0
	opts.on("-M INTEGER", "--minimun-coverage" "Set the minimun coverage to take a nucleotide position as mapped") do |m|
		options[:min_coverage] = m.to_i
	end

	options[:min_length] = 10
	opts.on("-L INTEGER", "--minimun-length" "Set the minimun length to take into account a enriched region. Default: 10") do |m|
		options[:min_length] = m.to_i
	end

	opts.on("-h", "--help", "Displays helps") do 
		puts opts
	end
end.parse!
#############################################################################################################################
## MAIN PROGRAM
#############################################################################################################################
control = load_file(options[:control], options[:min_coverage])
treatment = load_file(options[:treatment], options[:min_coverage])
enriched_chrs = compare_samples(control, treatment, options[:min_length])
save_output(enriched_chrs)