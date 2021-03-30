import subprocess
import os
import numpy as np

#########

n_reps = 1

#########

trees = ['three_taxon_all', 'three_taxon_medium', 'three_taxon_small', 'thirty_taxon_small', 'three_hundered_taxon_small']

directories = trees * n_reps
reps = np.repeat(range(1, n_reps+1), 5)

for x in range(0, len(reps)):
	directories[x] = directories[x] + '/' + str(reps[x])

grouped_directories = np.array_split(np.array(directories), n_reps)

number_of_datasets = [None]*n_reps
for x in range (0, n_reps):
	number_of_datasets[x] = [None] * 5

for x in directories:
	os.makedirs(os.path.join(x))

##########GENERATE_DIRECTORIES_FOR_MOLECULAR_DATA#####
	
dividers = ["directories <- c('"]  + ["', '"]*4

r_grouped_directories = np.array_split(np.array(directories), n_reps)	

final_data_groups = [None]*n_reps
for x in range (0, n_reps):
	final_data_groups[x] = [None] * 5

for x in range(0, n_reps):
	for y in range(0, 5):
		final_data_groups[x][y] = dividers[y] + r_grouped_directories[x][y]
	final_data_groups[x] = str(final_data_groups[x]).replace('", "', '').replace('["', '').replace('"]', "')")

####################
	
final_data_groups_rev = [None]*n_reps
for x in range (0, n_reps):
	final_data_groups_rev[x] = [None] * 5

for x in range(0, n_reps):
	for y in range(0, 5):
		final_data_groups_rev[x][y] = dividers[y] + r_grouped_directories[x][y]
	final_data_groups_rev[x] = str(final_data_groups_rev[x]).replace('", "', '').replace('["', '').replace('"]', "')").replace('c(', 'v(').replace("'", '"')
	
####################

scale_number = str([300, 30, 3, 30, 300]).replace('[', 'scale_number <- v(').replace(']', ')')
directory_denom = [2, 2, 2, 2, 24]

###########RUN_ANALYSES
	
for x in range(0, n_reps):
	R_script = open("r_shell.R", "r+")
	out_r_script = open("r_shell_analysis.R", "w+")
	Rscript_base = R_script.read() 
	out_r_script.write(final_data_groups[x] + Rscript_base)
	out_r_script.close()
	subprocess.call(['C:/Program Files/R/R-3.4.0/bin/Rscript', 'C:/Users/some3165/Documents/RATE_TIME_PROBLEM_FINAL/INDIVIDUAL_HETEROGENEITY/DIVERSIFICATION_RATE/ENTIRE_TREE/r_shell_analysis.R'])

	root_age_file = open('root_age.txt')
	for line in root_age_file:
		root_age = line.strip().split()	
	root_age = str(root_age).replace("['", "root_age <- ").replace("']", "")
	
	origin_age_file = open('origin_age.txt')
	for line in origin_age_file:
		origin_age = line.strip().split()	
	origin_age = str(origin_age).replace("['", "origin_age <- ").replace("']", "")
			
######SPECIFY_FOSSIL_DATA_FILES
		
	for y in range(4, 5):
		for z in range(1, (len(os.listdir(grouped_directories[x][y]))/9) +1):
			
###
			
			Oldest_fossil_amongst_all_fossils = open(str(grouped_directories[x][y]) + "/" + str(z) + "_r_output_oldest_fossil_ages.txt", "r+")
			Oldest_fossil_amongst_all_fossils_out = open(str(grouped_directories[x][y]) + "/" + str(z) + "_r_output_oldest_fossil_ages_out.txt", "w+")
						
###
			
			Fossil_ages = open(str(grouped_directories[x][y]) + "/" + str(z) + "_r_output_calibration_fossil_ages.txt", "r+")
			Fossil_ages_out = open(str(grouped_directories[x][y]) + "/" + str(z) + "_r_output_calibration_fossil_ages_out.txt", "w+")
						
###			
			
			Fossil_clades = open(str(grouped_directories[x][y]) + "/" + str(z) + "_taxon_sets.txt", "r+")
			Fossil_clades_out = open(str(grouped_directories[x][y]) + "/" + str(z) + "_taxon_sets_out.txt", "w+")
						
###
			
			Contains_calibrations = open(str(grouped_directories[x][y]) + "/" + str(z) + "_contains_calibrations.txt", "r+")
			Contains_calibrations_out = open(str(grouped_directories[x][y]) + "/" + str(z) + "_contains_calibrations_out.txt", "w+")
									
###
				
			Fossil_ages_base = Fossil_ages.read()
			numbers = range(0, 2000)
			for q in reversed(range(0, 2000)):
				Fossil_ages_base = Fossil_ages_base.replace("[[" + str(numbers[q]) + "]]", "")
			Fossil_ages_base = Fossil_ages_base.replace("\n\n\n[1] ", ", ").replace("\n[1] ", "fossil_ages <- v(").replace("\n", "") + ')'
			Fossil_ages_base = Fossil_ages_base.replace("\n\n", "")
			
			Oldest_fossil_amongst_all_fossils_base =  Oldest_fossil_amongst_all_fossils.read()
			for q in reversed(range(0, 2000)):
				Oldest_fossil_amongst_all_fossils_base = Oldest_fossil_amongst_all_fossils_base.replace("[" + str(numbers[q]) + "]", "[]")
			Oldest_fossil_amongst_all_fossils_base = Oldest_fossil_amongst_all_fossils_base.replace("[]", "fossil_ages_oldest <-")
			
			Fossil_clades_base = Fossil_clades.read()
			for q in reversed(range(0, 2000)):
				Fossil_clades_base = Fossil_clades_base.replace("[[" + str(numbers[q]) + "]]", "")
			for q in reversed(range(0, 2000)):	
				Fossil_clades_base = Fossil_clades_base.replace("[" + str(numbers[q]) + "]", "")
			Fossil_clades_base = Fossil_clades_base.replace('\n\n\n', '\ncalibration_taxa[++fossil_clade_counter] = clade(')
			Fossil_clades_base = Fossil_clades_base.replace(' ', '').replace('  ', '').replace('   ', '').replace('""', '", "').replace('\n"', ', "').replace('"\n', '")\n')
			Fossil_clades_base = 'fossil_clade_counter = 0\n' + Fossil_clades_base 
			Fossil_clades_base = Fossil_clades_base.replace('\n, ', '\ncalibration_taxa[++fossil_clade_counter] = clade(')
			Fossil_clades_base = Fossil_clades_base.replace('\n', '')
			Fossil_clades_base = Fossil_clades_base.replace('calibration_taxa', '\ncalibration_taxa')
			
			#
			
			Contains_calibrations_base = Contains_calibrations.read()
			Contains_calibrations_base = Contains_calibrations_base.replace('[1]', 'contains_calibrations <- ')
			
			###
			
			Fossil_ages_out.write(Fossil_ages_base)
			Oldest_fossil_amongst_all_fossils_out.write(Oldest_fossil_amongst_all_fossils_base)
			Fossil_clades_out.write(Fossil_clades_base)
			Contains_calibrations_out.write(Contains_calibrations_base)
			
			###
			
			Fossil_ages.close()
			Fossil_ages_out.close()
			Oldest_fossil_amongst_all_fossils.close()
			Oldest_fossil_amongst_all_fossils_out.close()
			Fossil_clades.close()
			Fossil_clades_out.close()
			Contains_calibrations.close()
			Contains_calibrations_out.close()
			
			###
			
			num_calibrated_nodes = -1
			with open(str(grouped_directories[x][y]) + "/" + str(z) + "_taxon_sets_out.txt", 'r') as Fossil_clades: 
				for line in Fossil_clades:
					num_calibrated_nodes += 1
						
			###
			
			Calibration_script_insertion = open("calibration_insertion.Rev", "r+")
			Calibration_script_insertion_out = open(str(grouped_directories[x][y]) + "/" + str(z) + "_calibration_insertion_out.Rev", "w+")
			
			Calibration_script_insertion_base = Calibration_script_insertion.read()
			Calibration_script_insertion_base = Calibration_script_insertion_base * num_calibrated_nodes
			Calibration_script_insertion_base = Calibration_script_insertion_base.replace(')tmrca', ')\n\ntmrca')
	
			Calibration_script_insertion_out.write('calibration_counter_two = 0\ncalibration_counter_three = 0\ncalibration_counter_four = 0\ncalibration_counter_five = 0\ncalibration_counter_six = 0\ncalibration_counter_seven = 0\ncalibration_counter_eight = 0\ncalibration_counter_nine = 0\ncalibration_counter_ten = 0\ncalibration_counter_eleven = 0\ncalibration_counter_twelve = 0\ncalibration_counter_thirteen = 0\ncalibration_counter_fourteen = 0\n\n' + Calibration_script_insertion_base)
			
			Calibration_script_insertion.close()
			Calibration_script_insertion_out.close()

			###
			
	for y in range(0, 5):
		number_of_datasets[x][y] = len(os.listdir(grouped_directories[x][y]))/directory_denom[y]
	dataset_number = str(number_of_datasets[x]).replace('[', 'dataset_number <- v(').replace(']', ')')
	
	Rev_script = open("rev_shell.Rev")
	out_rev_script = open("rev_shell_analysis.Rev", "w+")
	Revscript_base = Rev_script.read()
	out_rev_script.write(final_data_groups_rev[x] + '\n' + scale_number + '\n' + dataset_number + '\n' + root_age + '\n' + origin_age + '\n' + Revscript_base)
	Rev_script.close()
	out_rev_script.close()
	subprocess.call(['rb.exe', 'rev_shell_analysis.Rev'])