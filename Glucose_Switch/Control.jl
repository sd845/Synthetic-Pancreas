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
# Function: calculate_transcription_control_array
# Description: Calculate the transcriptional control array at time t
# Generated on: 2020-04-30T11:34:23.532
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_transcription_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control -
	control_array = zeros(6)

	# Alias the species -
	GFPssrA = x[1]
	GntR = x[2]
	RFPssrA = x[3]
	aplha_sigma28 = x[4]
	sigma_28 = x[5]
	sigma_70 = x[6]
	mRNA_GFPssrA = x[7]
	mRNA_GntR = x[8]
	mRNA_RFPssrA = x[9]
	mRNA_aplha_sigma28 = x[10]
	mRNA_sigma_28 = x[11]
	mRNA_sigma_70 = x[12]
	protein_GFPssrA = x[13]
	protein_GntR = x[14]
	protein_RFPssrA = x[15]
	protein_aplha_sigma28 = x[16]
	protein_sigma_28 = x[17]
	protein_sigma_70 = x[18]

	# Alias the binding parameters -
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	n_GFPssrA_sigma_70 = binding_parameter_dictionary["n_GFPssrA_sigma_70"]
	K_GFPssrA_sigma_70 = binding_parameter_dictionary["K_GFPssrA_sigma_70"]
	n_GFPssrA_GntR = binding_parameter_dictionary["n_GFPssrA_GntR"]
	K_GFPssrA_GntR = binding_parameter_dictionary["K_GFPssrA_GntR"]
	n_GntR_sigma_70 = binding_parameter_dictionary["n_GntR_sigma_70"]
	K_GntR_sigma_70 = binding_parameter_dictionary["K_GntR_sigma_70"]
	n_RFPssrA_sigma_28 = binding_parameter_dictionary["n_RFPssrA_sigma_28"]
	K_RFPssrA_sigma_28 = binding_parameter_dictionary["K_RFPssrA_sigma_28"]
	n_RFPssrA_aplha_sigma28 = binding_parameter_dictionary["n_RFPssrA_aplha_sigma28"]
	K_RFPssrA_aplha_sigma28 = binding_parameter_dictionary["K_RFPssrA_aplha_sigma28"]
	n_aplha_sigma28_sigma_70 = binding_parameter_dictionary["n_aplha_sigma28_sigma_70"]
	K_aplha_sigma28_sigma_70 = binding_parameter_dictionary["K_aplha_sigma28_sigma_70"]
	n_aplha_sigma28_GntR = binding_parameter_dictionary["n_aplha_sigma28_GntR"]
	K_aplha_sigma28_GntR = binding_parameter_dictionary["K_aplha_sigma28_GntR"]
	n_sigma_28_sigma_70 = binding_parameter_dictionary["n_sigma_28_sigma_70"]
	K_sigma_28_sigma_70 = binding_parameter_dictionary["K_sigma_28_sigma_70"]

	# Alias the control function parameters -
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_GFPssrA_RNAP = control_parameter_dictionary["W_GFPssrA_RNAP"]
	W_GFPssrA_sigma_70 = control_parameter_dictionary["W_GFPssrA_sigma_70"]
	W_GFPssrA_GntR = control_parameter_dictionary["W_GFPssrA_GntR"]
	W_GntR_RNAP = control_parameter_dictionary["W_GntR_RNAP"]
	W_GntR_sigma_70 = control_parameter_dictionary["W_GntR_sigma_70"]
	W_RFPssrA_RNAP = control_parameter_dictionary["W_RFPssrA_RNAP"]
	W_RFPssrA_sigma_28 = control_parameter_dictionary["W_RFPssrA_sigma_28"]
	W_RFPssrA_aplha_sigma28 = control_parameter_dictionary["W_RFPssrA_aplha_sigma28"]
	W_aplha_sigma28_RNAP = control_parameter_dictionary["W_aplha_sigma28_RNAP"]
	W_aplha_sigma28_sigma_70 = control_parameter_dictionary["W_aplha_sigma28_sigma_70"]
	W_aplha_sigma28_GntR = control_parameter_dictionary["W_aplha_sigma28_GntR"]
	W_sigma_28_RNAP = control_parameter_dictionary["W_sigma_28_RNAP"]
	W_sigma_28_sigma_70 = control_parameter_dictionary["W_sigma_28_sigma_70"]
	W_sigma_70_RNAP = control_parameter_dictionary["W_sigma_70_RNAP"]

# Gluconate binds to GntR like allolactose to LacI
# f_bound represents the fraction of bound Gluconate [Gluconate]^n/([Gluconate]^n+K^n)
# unbound GntR is represented by (1-f_bound)[GntR]
# Gluconate is passed as an parameter
gluconate_parameter_dictionary = data_dictionary["gluconate_parameter_dictionary"]
protein_gluconate = gluconate_parameter_dictionary["gluconate_concentration"]
n_gluconate_GntR = gluconate_parameter_dictionary["n_gluconate_GntR"]
K_gluconate_GntR = gluconate_parameter_dictionary["K_gluconate_GntR"]
# protein_GntR = x[14]

	f_bound = (protein_gluconate^(n_gluconate_GntR))/(protein_gluconate^(n_gluconate_GntR)+K_gluconate_GntR^(n_gluconate_GntR))
	protein_GntR = (1-f_bound)*protein_GntR

	# Transfer function target:GFPssrA actor:sigma_70
	actor_set_GFPssrA_sigma_70 = [
		protein_sigma_70
	]
	actor = prod(actor_set_GFPssrA_sigma_70)
	b_GFPssrA_sigma_70 = (actor^(n_GFPssrA_sigma_70))/(K_GFPssrA_sigma_70^(n_GFPssrA_sigma_70)+actor^(n_GFPssrA_sigma_70))

	# Transfer function target:GFPssrA actor:GntR
	actor_set_GFPssrA_GntR = [
		protein_GntR
	]
	actor = prod(actor_set_GFPssrA_GntR)
	b_GFPssrA_GntR = (actor^(n_GFPssrA_GntR))/(K_GFPssrA_GntR^(n_GFPssrA_GntR)+actor^(n_GFPssrA_GntR))

	# Control function for GFPssrA -
	control_array[1] = (W_GFPssrA_RNAP+W_GFPssrA_sigma_70*b_GFPssrA_sigma_70)/(1+W_GFPssrA_RNAP+W_GFPssrA_sigma_70*b_GFPssrA_sigma_70+W_GFPssrA_GntR*b_GFPssrA_GntR)

	# Transfer function target:GntR actor:sigma_70
	actor_set_GntR_sigma_70 = [
		protein_sigma_70
	]
	actor = prod(actor_set_GntR_sigma_70)
	b_GntR_sigma_70 = (actor^(n_GntR_sigma_70))/(K_GntR_sigma_70^(n_GntR_sigma_70)+actor^(n_GntR_sigma_70))

	# Control function for GntR -
	control_array[2] = (W_GntR_RNAP+W_GntR_sigma_70*b_GntR_sigma_70)/(1+W_GntR_RNAP+W_GntR_sigma_70*b_GntR_sigma_70)

	# Transfer function target:RFPssrA actor:sigma_28
	actor_set_RFPssrA_sigma_28 = [
		protein_sigma_28
	]
	actor = prod(actor_set_RFPssrA_sigma_28)
	b_RFPssrA_sigma_28 = (actor^(n_RFPssrA_sigma_28))/(K_RFPssrA_sigma_28^(n_RFPssrA_sigma_28)+actor^(n_RFPssrA_sigma_28))

	# Transfer function target:RFPssrA actor:aplha_sigma28
	actor_set_RFPssrA_aplha_sigma28 = [
		protein_aplha_sigma28
	]
	actor = prod(actor_set_RFPssrA_aplha_sigma28)
	b_RFPssrA_aplha_sigma28 = (actor^(n_RFPssrA_aplha_sigma28))/(K_RFPssrA_aplha_sigma28^(n_RFPssrA_aplha_sigma28)+actor^(n_RFPssrA_aplha_sigma28))

	# Control function for RFPssrA -
	control_array[3] = (W_RFPssrA_RNAP+W_RFPssrA_sigma_28*b_RFPssrA_sigma_28)/(1+W_RFPssrA_RNAP+W_RFPssrA_sigma_28*b_RFPssrA_sigma_28+W_RFPssrA_aplha_sigma28*b_RFPssrA_aplha_sigma28)

	# Transfer function target:aplha_sigma28 actor:sigma_70
	actor_set_aplha_sigma28_sigma_70 = [
		protein_sigma_70
	]
	actor = prod(actor_set_aplha_sigma28_sigma_70)
	b_aplha_sigma28_sigma_70 = (actor^(n_aplha_sigma28_sigma_70))/(K_aplha_sigma28_sigma_70^(n_aplha_sigma28_sigma_70)+actor^(n_aplha_sigma28_sigma_70))

	# Transfer function target:aplha_sigma28 actor:GntR

	actor_set_aplha_sigma28_GntR = [
		protein_GntR
	]
	actor = prod(actor_set_aplha_sigma28_GntR)
	b_aplha_sigma28_GntR = (actor^(n_aplha_sigma28_GntR))/(K_aplha_sigma28_GntR^(n_aplha_sigma28_GntR)+actor^(n_aplha_sigma28_GntR))

	# Control function for aplha_sigma28 -
	control_array[4] = (W_aplha_sigma28_RNAP+W_aplha_sigma28_sigma_70*b_aplha_sigma28_sigma_70)/(1+W_aplha_sigma28_RNAP+W_aplha_sigma28_sigma_70*b_aplha_sigma28_sigma_70+W_aplha_sigma28_GntR*b_aplha_sigma28_GntR)

	# Transfer function target:sigma_28 actor:sigma_70
	actor_set_sigma_28_sigma_70 = [
		protein_sigma_70
	]
	actor = prod(actor_set_sigma_28_sigma_70)
	b_sigma_28_sigma_70 = (actor^(n_sigma_28_sigma_70))/(K_sigma_28_sigma_70^(n_sigma_28_sigma_70)+actor^(n_sigma_28_sigma_70))

	# Control function for sigma_28 -
	control_array[5] = (W_sigma_28_RNAP+W_sigma_28_sigma_70*b_sigma_28_sigma_70)/(1+W_sigma_28_RNAP+W_sigma_28_sigma_70*b_sigma_28_sigma_70)

	# Control function for sigma_70 -
	control_array[6] = (W_sigma_70_RNAP)/(1+W_sigma_70_RNAP)

	# return -
	return control_array
end

#
# ----------------------------------------------------------------------------------- #
# Function: calculate_translation_control_array
# Description: Calculate the translation control array at time t
# Generated on: 2020-04-30T11:34:23.665
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Translation control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_translation_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control -
	control_array = ones(6)

	# return -
	return control_array
end
