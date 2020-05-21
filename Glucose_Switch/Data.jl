# ----------------------------------------------------------------------------------- #
# Copyright (c) 2020 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2020-04-30T11:35:41.902
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{String,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function build_data_dictionary(time_span::Tuple{Float64,Float64,Float64}, path_to_biophysical_constants_file::String = "./Default.json", host_type::Symbol = :bacteria)::Dict{String,Any}

	# load the biophysical_constants dictionary
	biophysical_constants_dictionary = build_biophysical_dictionary(path_to_biophysical_constants_file, host_type)

	# stoichiometric_matrix and dilution_matrix -
	stoichiometric_matrix = readdlm("./Network.dat")

	# number of states, and rates -
	(number_of_states,number_of_rates) = size(stoichiometric_matrix)

	# array of species types -
	species_symbol_type_array = [
		:gene	;	# 1	GFPssrA
		:gene	;	# 2	GntR
		:gene	;	# 3	RFPssrA
		:gene	;	# 4	aplha_sigma28
		:gene	;	# 5	sigma_28
		:gene	;	# 6	sigma_70
		:mrna	;	# 7	mRNA_GFPssrA
		:mrna	;	# 8	mRNA_GntR
		:mrna	;	# 9	mRNA_RFPssrA
		:mrna	;	# 10	mRNA_aplha_sigma28
		:mrna	;	# 11	mRNA_sigma_28
		:mrna	;	# 12	mRNA_sigma_70
		:protein	;	# 13	protein_GFPssrA
		:protein	;	# 14	protein_GntR
		:protein	;	# 15	protein_RFPssrA
		:protein	;	# 16	protein_aplha_sigma28
		:protein	;	# 17	protein_sigma_28
		:protein	;	# 18	protein_sigma_70
	]

	# we need to store the species symbol array for later -
	biophysical_constants_dictionary["species_symbol_type_array"] = species_symbol_type_array

	# array of gene lengths -
	gene_coding_length_array = [
		711.0	;	# 1	GFPssrA
		996.0	;	# 2	GntR
		714.0	;	# 3	RFPssrA
		720.0	;	# 4	aplha_sigma28
		720.0	;	# 5	sigma_28
		720.0	;	# 6	sigma_70
	]

	# array of mRNA coding lengths -
	mRNA_coding_length_array = [
		gene_coding_length_array[1]	;	# 7	1	mRNA_GFPssrA
		gene_coding_length_array[2]	;	# 8	2	mRNA_GntR
		gene_coding_length_array[3]	;	# 9	3	mRNA_RFPssrA
		gene_coding_length_array[4]	;	# 10	4	mRNA_aplha_sigma28
		gene_coding_length_array[5]	;	# 11	5	mRNA_sigma_28
		gene_coding_length_array[6]	;	# 12	6	mRNA_sigma_70
	]

	# array of mRNA coding lengths -
	protein_coding_length_array = [
		round((0.33)*mRNA_coding_length_array[1])	;	# 13	1	protein_GFPssrA
		round((0.33)*mRNA_coding_length_array[2])	;	# 14	2	protein_GntR
		round((0.33)*mRNA_coding_length_array[3])	;	# 15	3	protein_RFPssrA
		round((0.33)*mRNA_coding_length_array[4])	;	# 16	4	protein_aplha_sigma28
		round((0.33)*mRNA_coding_length_array[5])	;	# 17	5	protein_sigma_28
		round((0.33)*mRNA_coding_length_array[6])	;	# 18	6	protein_sigma_70
	]

	# array of gene concentrations -
	gene_abundance_array = [
		5.0	;	# (nM) 1	GFPssrA
		5.0	;	# (nM) 2	GntR
		5.0	;	# (nM) 3	RFPssrA
		5.0	;	# (nM) 4	aplha_sigma28
		5.0	;	# (nM) 5	sigma_28
		0.0	;	# (nM) 6	sigma_70
	]

	# initial condition array -
	initial_condition_array = [
		gene_abundance_array[1]	;	# 1	GFPssrA
		gene_abundance_array[2]	;	# 2	GntR
		gene_abundance_array[3]	;	# 3	RFPssrA
		gene_abundance_array[4]	;	# 4	aplha_sigma28
		gene_abundance_array[5]	;	# 5	sigma_28
		gene_abundance_array[6]	;	# 6	sigma_70
		0.0	;	# 7	mRNA_GFPssrA
		0.0	;	# 8	mRNA_GntR
		0.0	;	# 9	mRNA_RFPssrA
		0.0	;	# 10	mRNA_aplha_sigma28
		0.0	;	# 11	mRNA_sigma_28
		0.0	;	# 12	mRNA_sigma_70
		0.0	;	# 13	protein_GFPssrA
		0.0	;	# 14	protein_GntR
		0.0	;	# 15	protein_RFPssrA
		0.0	;	# 16	protein_aplha_sigma28
		0.0	;	# 17	protein_sigma_28
		0.035	;	# 18	protein_sigma_70
	]

	binding_parameter_dictionary = Dict{String,Float64}()
	binding_parameter_dictionary["n_GFPssrA_sigma_70"] = 0.8
	binding_parameter_dictionary["K_GFPssrA_sigma_70"] = 38.14
	binding_parameter_dictionary["n_GFPssrA_GntR"] = 1.0
	binding_parameter_dictionary["K_GFPssrA_GntR"] = 30.0
	binding_parameter_dictionary["n_GntR_sigma_70"] = 1.0
	binding_parameter_dictionary["K_GntR_sigma_70"] = 30.0
	binding_parameter_dictionary["n_RFPssrA_sigma_28"] = 1.0
	binding_parameter_dictionary["K_RFPssrA_sigma_28"] = 30.0
	binding_parameter_dictionary["n_RFPssrA_aplha_sigma28"] = 1.7
	binding_parameter_dictionary["K_RFPssrA_aplha_sigma28"] = 60.0
	binding_parameter_dictionary["n_aplha_sigma28_sigma_70"] = 0.52
	binding_parameter_dictionary["K_aplha_sigma28_sigma_70"] = 46.8
	binding_parameter_dictionary["n_aplha_sigma28_GntR"] = 1.0
	binding_parameter_dictionary["K_aplha_sigma28_GntR"] = 30.0
	binding_parameter_dictionary["n_sigma_28_sigma_70"] = 0.52
	binding_parameter_dictionary["K_sigma_28_sigma_70"] = 46.8

	# Alias the control function parameters -
	control_parameter_dictionary = Dict{String,Float64}()
	control_parameter_dictionary["W_GFPssrA_RNAP"] = 0.000014
	control_parameter_dictionary["W_GFPssrA_sigma_70"] = 10.0
	control_parameter_dictionary["W_GFPssrA_GntR"] = 1.0
	control_parameter_dictionary["W_GntR_RNAP"] = 0.001
	control_parameter_dictionary["W_GntR_sigma_70"] = 1.0
	control_parameter_dictionary["W_RFPssrA_RNAP"] = 0.001
	control_parameter_dictionary["W_RFPssrA_sigma_28"] = 1.0
	control_parameter_dictionary["W_RFPssrA_aplha_sigma28"] = 100#1.0
	control_parameter_dictionary["W_aplha_sigma28_RNAP"] = 0.001
	control_parameter_dictionary["W_aplha_sigma28_sigma_70"] = 1.0
	control_parameter_dictionary["W_aplha_sigma28_GntR"] = 1.0
	control_parameter_dictionary["W_sigma_28_RNAP"] = 0.00014
	control_parameter_dictionary["W_sigma_28_sigma_70"] = 100.0
	control_parameter_dictionary["W_sigma_70_RNAP"] = 0.0

	# degradation modifiers -
	degradation_modifier_array = [
		0.0	;	# 1	GFPssrA
		0.0	;	# 2	GntR
		0.0	;	# 3	RFPssrA
		0.0	;	# 4	aplha_sigma28
		0.0	;	# 5	sigma_28
		0.0	;	# 6	sigma_70
		0.05	;	# 7	mRNA_GFPssrA
		1.0	;	# 8	mRNA_GntR
		1.0	;	# 9	mRNA_RFPssrA
		0.05	;	# 10	mRNA_aplha_sigma28
		0.05	;	# 11	mRNA_sigma_28
		1.0	;	# 12	mRNA_sigma_70
		28.0	;	# 13	protein_GFPssrA
		1.0	;	# 14	protein_GntR
		30.0	;	# 15	protein_RFPssrA
		1.0	;	# 16	protein_aplha_sigma28
		1.0	;	# 17	protein_sigma_28
		1.0	;	# 18	protein_sigma_70
	]


	# time constant modifiers -
	time_constant_modifier_array = [
		0.0	;	# 1	GFPssrA
		0.0	;	# 2	GntR
		0.0	;	# 3	RFPssrA
		0.0	;	# 4	aplha_sigma28
		0.0	;	# 5	sigma_28
		0.0	;	# 6	sigma_70
		1.0	;	# 7	mRNA_GFPssrA
		1.0	;	# 8	mRNA_GntR
		1.0	;	# 9	mRNA_RFPssrA
		1.0	;	# 10	mRNA_aplha_sigma28
		1.0	;	# 11	mRNA_sigma_28
		1.0	;	# 12	mRNA_sigma_70
		1.0	;	# 13	protein_GFPssrA
		1.0	;	# 14	protein_GntR
		1.0	;	# 15	protein_RFPssrA
		1.0	;	# 16	protein_aplha_sigma28
		1.0	;	# 17	protein_sigma_28
		1.0	;	# 18	protein_sigma_70
	]

	# Dilution degrdation matrix -
	dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)

	# Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array, host_type)

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)

	# D- Gluconate parameter values
	gluconate_parameter_dictionary = Dict{String,Float64}()
	gluconate_parameter_dictionary["gluconate_concentration"] = 0.0
	gluconate_parameter_dictionary["n_gluconate_GntR"] = 2.45
	gluconate_parameter_dictionary["K_gluconate_GntR"] = 5.3 #muM

	# Parameter name index array -
	parameter_name_mapping_array = [
		"n_GFPssrA_sigma_70"	;	# 1
		"K_GFPssrA_sigma_70"	;	# 2
		"n_GFPssrA_GntR"	;	# 3
		"K_GFPssrA_GntR"	;	# 4
		"n_GntR_sigma_70"	;	# 5
		"K_GntR_sigma_70"	;	# 6
		"n_RFPssrA_sigma_28"	;	# 7
		"K_RFPssrA_sigma_28"	;	# 8
		"n_RFPssrA_aplha_sigma28"	;	# 9
		"K_RFPssrA_aplha_sigma28"	;	# 10
		"n_aplha_sigma28_sigma_70"	;	# 11
		"K_aplha_sigma28_sigma_70"	;	# 12
		"n_aplha_sigma28_GntR"	;	# 13
		"K_aplha_sigma28_GntR"	;	# 14
		"n_sigma_28_sigma_70"	;	# 15
		"K_sigma_28_sigma_70"	;	# 16
		"W_GFPssrA_RNAP"	;	# 17
		"W_GFPssrA_sigma_70"	;	# 18
		"W_GFPssrA_GntR"	;	# 19
		"W_GntR_RNAP"	;	# 20
		"W_GntR_sigma_70"	;	# 21
		"W_RFPssrA_RNAP"	;	# 22
		"W_RFPssrA_sigma_28"	;	# 23
		"W_RFPssrA_aplha_sigma28"	;	# 24
		"W_aplha_sigma28_RNAP"	;	# 25
		"W_aplha_sigma28_sigma_70"	;	# 26
		"W_aplha_sigma28_GntR"	;	# 27
		"W_sigma_28_RNAP"	;	# 28
		"W_sigma_28_sigma_70"	;	# 29
		"W_sigma_70_RNAP"	;	# 30
		"rnapII_concentration"	;	# 31
		"ribosome_concentration"	;	# 32
		"degradation_constant_mRNA"	;	# 33
		"degradation_constant_protein"	;	# 34
		"kcat_transcription"	;	# 35
		"kcat_translation"	;	# 36
		"maximum_specific_growth_rate"	;	# 37
		"saturation_constant_transcription"	;	# 38
		"saturation_constant_translation"	;	# 39
	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{String,Any}()
	data_dictionary["number_of_states"] = number_of_states
	data_dictionary["species_symbol_type_array"] = species_symbol_type_array
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["gene_coding_length_array"] = gene_coding_length_array
	data_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	data_dictionary["protein_coding_length_array"] = protein_coding_length_array
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix
	data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
	data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
	data_dictionary["parameter_name_mapping_array"] = parameter_name_mapping_array
	data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array
	data_dictionary["translation_parameter_array"] = translation_parameter_array
	data_dictionary["degradation_modifier_array"] = degradation_modifier_array
	data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array
	data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary
	data_dictionary["gluconate_parameter_dictionary"] = gluconate_parameter_dictionary
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	data_dictionary["R"] = 8.314 			# J mol^-1 K^-1
	data_dictionary["T_K"] = 273.15 + 29.0 	# K

	return data_dictionary
end
