#! /usr/bin/env ruby

require 'optparse'


##############################################################################################################################
## METHODS
##############################################################################################################################
def load_and_clean_file(files_array, ind_seq_size)  #need to create a empty hash (hash_of_sequences) 
	all_fastas = {}
	indeterminated_seq = "N" * ind_seq_size
	files_array.each do |file|	
		sequences = {}
		name = nil
		seq = ''
		File.open(file).each do |line|  #opens the file and builds the hash "sequences" with the fasta sequences
			line.chomp!
			if line =~ /^>/
				seq = '' if !indeterminated_seq.empty? && seq.include?(indeterminated_seq) 
				sequences[name.split(/\s/).first] = seq if !name.nil? && !seq.empty?
				seq = ''
				name = line.gsub('>', '')
			else
				seq << line
			end 
		end
		sequences[name.split(/\s/).first] = seq if !name.nil? && !seq.empty?
		all_fastas[file] = sequences
	end
	return all_fastas
end

def discard_seqs_from_list(all_fastas, discarded_list_filename)
	sequences_to_discard = load_list(discarded_list_filename)
	all_fastas.each do |filename, sequences|
		sequences_to_discard.each do |sequence|
			sequences.delete(sequence)
		end
	end
end

def select_seqs_from_list(all_fastas, acepted_list_filename)
	sequences_to_acept = load_list(acepted_list_filename)
	all_fastas.each do |filename, sequences|
		sequences.reject! { |name, seq| !sequences_to_acept.include?(name) }
	end
end

def load_list(filename)
	list = []
	File.open(filename).each do |line|
		line.chomp!
		list << line
	end
	return list
end

def filter_larger_seqs(all_fastas, filter, output_file) #all_fastas is a hash of hashes 
	filtered_fastas = {}
	discarded_fastas = {}
	all_fastas.each do |filename, fastas| 
		filtered_seqs = {}
		discarded_seqs = {}
		fastas.each do |id, seq|
			if seq.size >= filter
				discarded_seqs[id] = seq
			elsif seq.size < filter
				filtered_seqs[id] = seq
			end
		end
		discarded_fastas[filename] = discarded_seqs if !discarded_seqs.empty?
		filtered_fastas[filename] = filtered_seqs if !filtered_seqs.empty?
	end
	save_output_in_a_file("larger_seqs_#{output_file}", discarded_fastas)
	return filtered_fastas
end

def filter_smaller_seqs(all_fastas, filter, output_file) #all_fastas is a hash of hashes 
	filtered_fastas = {}
	discarded_fastas = {}
	all_fastas.each do |filename, fastas| 
		filtered_seqs = {}
		discarded_seqs = {}
		fastas.each do |id, seq|
			if seq.size <= filter
				discarded_seqs[id] = seq
			elsif seq.size > filter
				filtered_seqs[id] = seq
			end
		end
		discarded_fastas[filename] = discarded_seqs if !discarded_seqs.empty?
		filtered_fastas[filename] = filtered_seqs if !filtered_seqs.empty?
	end
	save_output_in_a_file("smaller_seqs_#{output_file}", discarded_fastas)
	return filtered_fastas
end

def get_per_sequence_statistics(all_fastas)
	all_fastas.each do |filename, fastas|
		puts filename
		fastas.each do |id, seq|
			gc_content = nil
			type = nil
			if seq.downcase.include?("u")
				gc_content = get_per_GC(seq)
				type = 'rna'
			elsif seq.downcase =~ /[^agtcn]/
				type = 'prot'
			else
				gc_content = get_per_GC(seq)
				type = 'dna'
			end
			puts "#{id}\t#{type}\t#{seq.size}\t#{gc_content}"
		end
	end
	puts "\n"
end

def get_per_GC(x)
	sum_gc = x.downcase.count('g') + x.downcase.count('c')
	return (sum_gc.fdiv(x.size) * 100).round(5)
end

def load_motifs(motif_string) 
	motifs = []
	if File.exist?(motif_string) #input can be a file or a string
		File.open(motif_string).each do |line| 
			line.chomp!
			motifs << line.upcase if !line.empty?
		end
	else
		motif_string.split(',').each do |motif|
			motifs << motif.upcase
		end
	end
	return motifs   
	
end

def find_motifs(all_files, motifs) #sequences and fastas_motifs must be a hash of hashes and motifs must be an array
	fastas_motifs = {}
	if !motifs.empty?
		all_files.each do |file, sequences|
			sequences.each do | id, seq|
				seq_match = {}
				motifs.each do |motif|
					position = 0
					matches = []
					
					while !position.nil?
						position = seq.downcase.index(motif.downcase, position) # return based 0 coordinates
						if !position.nil?
							position += 1 # save based 1 coordinates
							matches << position
						end
					end
					seq_match[motif] = matches if !matches.empty?
				end
				fastas_motifs[id] = seq_match if !seq_match.empty?
			end
		end
	end
	return fastas_motifs
end

def get_n50(all_fastas)	
	seq_lenghts = []	#all_fastas -> hash of hashes
	all_fastas.each do |filename, fastas|
		fastas.each do |id, seq|
			seq_lenghts << seq.size
		end
	end
	seq_lenghts.sort!
	n50 = 0
	if seq_lenghts.count.odd?
		median = seq_lenghts.count / 2
		puts "N50: #{seq_lenghts[median]}"
	else
		median = (seq_lenghts.count + 1) / 2
		puts "N50: #{seq_lenghts[median]}"
	end
end

def get_max_and_min_GC_content(all_fastas)                                             
	all_fastas.each do |file_name, sequences|
		max = 0
		min = 100
		name_max = nil
		name_min = nil
		sequences.each do |name, seq|
			percentage_GC = get_per_GC(seq)
			if percentage_GC > max 
				max = percentage_GC
				name_max = name
			elsif percentage_GC < min
				min = percentage_GC
				name_min = name
			end
		end
		puts "File: #{file_name}\nThe sequence with the highest GC content --> \n#{name_max}\n#{max}\nThe sequence with the lowest GC content -->\n#{name_min}\n#{min}"
	end
end

def load_fragments_table(files)
	fragments = []
	files.each do |file|
		File.open(file).each do |line|
			line.chomp!
			fragments << line.split("\t") if !line.empty?
		end
	end
	return fragments
end

def find_fragments(all_fastas, fragments, tail_length) #fragments is an array of arrays
	output_fasta = {}
	all_fastas.each do |filename, fastas|
		output_seq = {}
		fastas.each do |id, seq|
			fragment_number = 0
			fragments.each do |name, start, stop|
				if id == (name)
					fragment_number += 1
					starting = start.to_i - tail_length - 1
					starting = 0 if starting < 0
					ending = stop.to_i + tail_length - 1
					ending = seq.size - 1 if ending > seq.size
					seq_name = "#{name}:#{starting + 1}:#{ending + 1}_seq_num_fasta_editor_#{fragment_number}"
					output_seq[seq_name] = seq[starting..ending]
				end
			end
		end
		output_fasta[filename] = output_seq  if !output_seq.empty?
	end
	return output_fasta
end

def generate_new_id(main_operating_hash_of_sequences)
	new_main_operating_hash_of_sequences = {}
	seq_number = 1
	main_operating_hash_of_sequences.each do |filename, fastas| 
		new_fastas = {}
		fastas.each do |id, seq|
			new_seq_name = "seq_#{seq_number}"
			new_fastas[new_seq_name] = seq
			seq_number += 1
			puts "#{id}\t#{new_seq_name}"
		end
		new_main_operating_hash_of_sequences[filename] = new_fastas
	end
	puts "\n"
	return new_main_operating_hash_of_sequences
end



## output methods 
#-----------------------------------------------------------------------------------------------------------


def format_all_fastas(all_fastas) #all_fastas is a hash of hashes 
	output_fasta = {}
	all_fastas.each do |filename, fastas| 
		output_seq = {}
		fastas.each do |id, seq|
			output_seq[id] = seq if !output_seq[id]
		end
		output_fasta[filename] = output_seq if !output_seq.empty?
	end
	return output_fasta
end

def format_output_with_sequences_with_motifs(all_fastas, fastas_motifs) 
	output_fasta = {}
	all_fastas.each do |filename, fastas|
		output_seq = {}
		fastas.each do |id, seq|
			fastas_motifs.each do |name_seq, matches|
				if id == name_seq || !output_seq[id]
					output_seq[id] = seq
				end
			end
		end
		output_fasta[filename] = output_seq if !output_seq.empty?
	end
	return output_fasta
end

def format_output_with_sequences_with_motifs_header(all_fastas, fastas_motifs)  ##############3333
	output_fasta = {}
	all_fastas.each do |filename, fastas|
		output_seq = {}
		fastas.each do |id, seq|
			fastas_motifs.each do |name_seq, matches|
				key_string = "#{id} " if id == name_seq
				matches.each do |motif, places|
					key_string.insert(-1, "#{motif}-#{places.join(',')}_") if key_string 
				end
				output_seq[key_string] = seq if !output_seq[key_string] 
			end
		end
		output_fasta[filename] = output_seq if !output_seq.empty?
	end
	return output_fasta
end

def save_output_in_a_file(filename, output_fasta)
	File.open(filename, 'w') do |file|
		output_fasta.each do |filename, fastas|
			fastas.each do |id, seq|
				file << ">#{id.split('_seq_num_fasta_editor_')[0]}\n#{seq}\n" if !id.nil?
			end
		end
	end
end



#############################################################################################################################
## INPUT PARSING
##################################################################################################################################3

options = {}

OptionParser.new do  |opts|
	options[:input] = []
	opts.on("-i FILE1,FILE2,FILE3", "--file", "Open file") do |i|
		options[:input] = i.split(',')		
	end

	options[:clean] = 0
	opts.on("-C INTEGER", "--clean", "Clean fasta removing sequences with gaps of FLOAT indeterminations") do |integer|
		options[:clean] = integer.to_i
	end

	options[:fragments] = nil
	opts.on("-f FILE1,FILE2,FILE3", "--fragments", "Find fragments in fasta sequences") do |f|
		options[:fragments] = f.split(",")
	end
	options[:tail] = 0
	opts.on("-T INTERGER", "--tail", "Length of the tails of the fragments. Require '-f'. [DEFAULT = 0]") do |integer|
		options[:tail] = integer.to_i
	end

	options[:filter_larger_seqs] = ''
	opts.on("-L INTEGER", "--filter-larger", "Discard sequences larger than INTEGER, and save it in 'FILENAME_larger_seqs.fasta'. Then the script works only whith the smaller sequences") do |integer|
		options[:filter_larger_seqs] = integer
	end

	options[:filter_smaller_seqs] = ''
	opts.on("-S INTEGER", "--filter-smaller", "Discard sequences smaller than INTEGER, and save it in 'FILENAME_smaller_seqs.fasta, Then the script works only whith the larger sequences' ") do |integer|
		options[:filter_smaller_seqs] = integer
	end

	options[:sequences_to_discard] = nil
	opts.on("-d FILE", "--discard-sequences", "Discard sequences from list") do |filename|
		options[:sequences_to_discard] = filename
	end
	options[:sequences_to_accept] = nil
	opts.on("-l FILE", "--list-sequences", "Only work with sequences on list") do |filename|
		options[:sequences_to_accept] = filename
	end

	options[:rename] = false
	opts.on("-r", "--rename", "Rename sequences") do
		options[:rename] = true
	end

	options[:sequence_statistics] = false
	opts.on("-s", "--sequence_statistics", "Summary table with statistics per each sequence") do 
		options[:sequence_statistics] = true
	end

	options[:content] = false
	opts.on("-c", "--content", "Return GC content of each sequence") do
		options[:content] = true
	end

	opts.on("-m MOTIF1,MOTIF2,MOTIF3", "--motif", "Find a motif") do |m|
		options[:motif] = m
	end

	options[:n50] = false
	opts.on("-n", "--n50", "Return N50 value") do 
		options[:n50] = true
	end

	options[:translate] = false
	opts.on("-t", "--translate", "Translate the sequences") do 
		options[:translate] = true
	end


	options[:create] = '' 
	opts.on("-c STRING", "--create", "Create an output: 'a' for export all fastas of all files, 'm'  export fastas with the matches of the motifs(require -m)
						'mh'  export fastas with the matches of the motifs showing them in the header(require -m)") do |c|
		options[:create] << c
	end
	
	options[:output] = ''
	opts.on("-o FILENAME", "--save", "save in a file with other name (require -c)") do |o|
		options[:output] = "#{o}"
	end

	opts.on("-h", "--help", "Displays helps") do 
		puts opts
	end

end.parse!

#############################################################################################################################
## MAIN PROGRAM
##################################################################################################################################3

abort("ERROR: The specified files not exist") if options[:input].empty?

all_fastas = load_and_clean_file(options[:input], options[:clean])
all_fastas = discard_seqs_from_list(all_fastas, options[:sequences_to_discard]) if !options[:sequences_to_discard].nil?
all_fastas = select_seqs_from_list(all_fastas, options[:sequences_to_accept]) if !options[:sequences_to_accept].nil?
all_fastas = filter_larger_seqs(all_fastas, options[:filter_larger_seqs].to_i, options[:output]) if !options[:filter_larger_seqs].empty?

all_fastas = filter_smaller_seqs(all_fastas, options[:filter_smaller_seqs].to_i, options[:output]) if !options[:filter_smaller_seqs].empty?

if !options[:fragments].nil?
	fragments = load_fragments_table(options[:fragments])
	sequences_to_process =  find_fragments(all_fastas, fragments, options[:tail])
else
	sequences_to_process = all_fastas
end

sequences_to_process = generate_new_id(sequences_to_process) if options[:rename]

get_per_sequence_statistics(sequences_to_process) if options[:sequence_statistics]
get_max_and_min_GC_content(sequences_to_process) if options[:content]
get_n50(sequences_to_process) if options[:n50]

fastas_motifs = {}
if options[:motif]
	motifs = load_motifs(options[:motif])
	fastas_motifs = find_motifs(sequences_to_process, motifs)
	fastas_motifs.each do |id, motifs|
		motifs.each do |motif, matches|
			puts "#{id}\t#{motif}\t#{matches.join(',')}"
		end
	end	
end


if options[:create].include?('a')
	output_fasta = sequences_to_process 
elsif options[:create].include?('m')
	if options[:create].include?('h')
		output_fasta = format_output_with_sequences_with_motifs_header(sequences_to_process, fastas_motifs)
	else
		output_fasta = format_output_with_sequences_with_motifs(sequences_to_process, fastas_motifs) 
	end
end

save_output_in_a_file(options[:output], output_fasta) if !options[:output].empty?

	










		
 
