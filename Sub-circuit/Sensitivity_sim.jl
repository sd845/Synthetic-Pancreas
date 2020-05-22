# include the includes -
include("Include.jl")

# mean center -
function mean_center_array(results_array::Array{Float64,2})::Array{Float64,2}

    # get the size -
    (NR,NC) = size(results_array)
    scaled_array = zeros(NR,NC)

    for col_index = 1:NC

        data_col = results_array[:,col_index]
        mu_value = mean(data_col)
        std_value = std(data_col)

        for row_index = 1:NR
            scaled_array[row_index,col_index] = (data_col[row_index] - mu_value)/(std_value)
        end
    end

    return scaled_array
end

# computes the model performance -
function model_performance(parameter_guess_array,index)

    # what is the host_type?
    host_type = :cell_free

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 960.0
    time_step_size = 0.01

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"
    #path_to_data_dir = "$(pwd())/data"

    # Load the data dictionary (uses the default biophysical_constants file)
    model_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # update the paramaters in the model data dictionary -
    # for now - lets only search over dG's -
    R = model_data_dictionary["R"]
    T_K = model_data_dictionary["T_K"]

    # what is the size of the parameter_guess_array?
    number_of_parameters = length(parameter_guess_array)

    # Phase 1: parameter update =========================================================================== #
    # update the paramaters in the model data dictionary -
    # for now - lets only search over dG's -
    R = model_data_dictionary["R"]
    T_K = model_data_dictionary["T_K"]

    # compute W -
    tmp_W_array = Float64[]
    for index = 1:5
        parameter_guess = parameter_guess_array[index]
        value = exp(-1*parameter_guess/(R*T_K))
        push!(tmp_W_array,value)
    end

    # update the control W's -
	control_parameter_dictionary = model_data_dictionary["control_parameter_dictionary"]
	control_parameter_dictionary["W_GntR_RNAP"] = tmp_W_array[1]
	control_parameter_dictionary["W_GntR_sigma_70"] = tmp_W_array[2]
	control_parameter_dictionary["W_deGFP_RNAP"] = tmp_W_array[3]
	control_parameter_dictionary["W_deGFP_sigma_70"] = tmp_W_array[4]
	control_parameter_dictionary["W_deGFP_GntR"] = tmp_W_array[5]
	#control_parameter_dictionary["W_sigma_70_RNAP"] = tmp_W_array[6]
	model_data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary

    binding_parameter_dictionary = model_data_dictionary["binding_parameter_dictionary"]
	binding_parameter_dictionary["n_GntR_sigma_70"] = parameter_guess_array[6]
	binding_parameter_dictionary["K_GntR_sigma_70"] = parameter_guess_array[7]
	binding_parameter_dictionary["n_deGFP_sigma_70"] = parameter_guess_array[8]
	binding_parameter_dictionary["K_deGFP_sigma_70"] = parameter_guess_array[9]
	binding_parameter_dictionary["n_deGFP_GntR"] = parameter_guess_array[10]
	binding_parameter_dictionary["K_deGFP_GntR"] = parameter_guess_array[11]
    model_data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

    # time constant modifier -
	time_constant_modifier_array = [
			0.0					;	# 1	GntR
			0.0					;	# 2	deGFP
			0.0					;	# 3	sigma_70
			parameter_guess_array[12]	;	# 4	mRNA_GntR #unknown
			parameter_guess_array[13] 	;	# 5	mRNA_deGFP
			1.0	;							# 6	mRNA_sigma_70
			parameter_guess_array[14]	;	# 7	protein_GntR #unknown
			parameter_guess_array[15]	;	# 8	protein_deGFP
			1.0	;							# 9	protein_sigma_70
	]
    model_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

    # setup degradation_modifier_array -
	degradation_modifier_array = [
			0.0	;	# 1	GntR #unknown
			0.0	;	# 2	deGFP
			0.0	;	# 3	sigma_70
			parameter_guess_array[16]	;	# 4	mRNA_GntR #unknown
			parameter_guess_array[17]	;	# 5	mRNA_deGFP
			0.0							;	# 6	mRNA_sigma_70
			parameter_guess_array[18]	;	# 7	protein_GntR #unknown
			parameter_guess_array[19]	;	# 8	protein_deGFP
			parameter_guess_array[20]	;	# 9	protein_sigma_70
	]
    model_data_dictionary["degradation_modifier_array"] = degradation_modifier_array

    # update the translation time -
    #model_data_dictionary["half_life_translation_capacity"] = parameter_guess_array[32]

    # lastly, update KL -
    biophysical_constants_dictionary = model_data_dictionary["biophysical_constants_dictionary"]
    #biophysical_constants_dictionary["translation_saturation_constant"] = parameter_guess_array[33]
    model_data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

    # grab defaults -
    species_symbol_type_array = model_data_dictionary["species_symbol_type_array"]
    protein_coding_length_array = model_data_dictionary["protein_coding_length_array"]
    gene_coding_length_array = model_data_dictionary["gene_coding_length_array"]
    time_constant_modifier_array = model_data_dictionary["time_constant_modifier_array"]
    initial_condition_array = model_data_dictionary["initial_condition_array"]

    # # get gene IC -
    idx_gene = findall(x->x==:gene,species_symbol_type_array)
    gene_abundance_array = initial_condition_array[idx_gene]

    # Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
    model_data_dictionary["translation_parameter_array"] = translation_parameter_array

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
    model_data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

    # Dilution degrdation matrix -
    dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)
    model_data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix
    # ===================================================================================================== #

    # Phase 2:  solve model equations ===================================================================== #
    # solve the balance equations -
    (TSIM,XSIM) = SolveBalances(time_start,time_stop,time_step_size,model_data_dictionary)
	#display(XSIM)
	#display(min(TSIM))
	#display(max(TSIM))
    # ===================================================================================================== #

    # Phase 3: compute the model performance metrics ====================================================== #
    p_GFP_AUC = integrate(TSIM,XSIM[:,index],SimpsonEven())
    # ===================================================================================================== #

    # return the performance_array -
    return p_GFP_AUC

end

function main(index)

    # setup the sensitivity function -
	#P = zeros(51)
    SF(P) = model_performance(P,index)

	parameter_bounds_array = Array{Tuple,1}()
	size(parameter_bounds_array)
    # setup ranges -
    #sample_bounds_array = Array{Tuple,1}()
	sample_bounds_array = [

	        # dG's -
	        39995.69565886529 68282.2841495398      ;   # 1 W_GntR_RNAP
	        -50000.0 25000.0  						;   # 2 W_GntR_sigma_70
			39995.69565886529 68282.2841495398  	;   # 3 W_deGFP_RNAP
			-50000.0 25000.0  						;	# 4 W_deGFP_sigma_70
			-50000.0 25000.0   			  			;   # 5 W_deGFP_GntR
			# W_sigma_70_RNAP = 0.0

	        # binding parameters -
			0.5 10.0 				   	; #6  n_GntR_sigma_70
			0.001 100.0					; #7  K_GntR_sigma_70
			0.5037004212911795 2.2  	; #8  n_deGFP_sigma_70
			25.539812443045424 100.0 	; #9  K_deGFP_sigma_70
			0.5 10.0         			; #10 n_deGFP_GntR
			0.001 100.0         		; #11 K_deGFP_GntR

	        # time constants -
			0.0001 100.0							;	# 12	mRNA_GntR #unknown
			0.0039273045649573235 7.17667775643900 	;	# 13	mRNA_deGFP
			0.0001 100.0							;	# 14	protein_GntR #unknown
			0.35035678021027006 10.203484322453043  ;	# 15	protein_deGFP

	        # degradation mods -
			0.1 10.0	;	# 16	mRNA_GntR #unknown
			0.1 10.0	;	# 17	mRNA_deGFP
			0.1 10.0	;	# 18	protein_GntR #unknown
			40.0 60.0	;	# 19	protein_deGFP
			0.1 10.0	;	# 20	protein_sigma_70

	    ];

    #@show index

	number_of_parameters = size(sample_bounds_array)[1]
			for parameter_index = 1:(number_of_parameters)
				display(parameter_index)
				# get row of parameters -
				lb = sample_bounds_array[parameter_index,1]
				ub = sample_bounds_array[parameter_index,2]

				#display(lb)
				#display(ub)
				# create the tuple -
				tmp_tuple = (lb,ub)
				#display(tmp_tuple)
				# cache -
				push!(parameter_bounds_array,tmp_tuple)
				#parameter_bounds_array[parameter_index] = tmp_tuple
				display(parameter_bounds_array)
			end
#size(parameter_bounds_array)
    # do the global sensitivity analysis -
    sensitivity_results = gsa(SF,Morris(total_num_trajectory=10000,num_trajectory=1000),parameter_bounds_array)

    # return -
    return sensitivity_results
end

# setup paths -
#path_to_ensemble_file = "$(pwd())/Ensemble-c1-restriction-T20.dat"

# compute a sensitivity array for the AUC of each species -
species_index_array = [4 5 7 8]
number_of_species = length(species_index_array)
number_of_parameters = 20
results_array = zeros(number_of_parameters,1)
for species_index in species_index_array

    global results_array

    # conduct senstivity analysis -
    sensitivity_results = main(species_index)

    # get the μ and σ^2
    mu = sensitivity_results.means
    var = sensitivity_results.variances

    #@show mu, var
    results_array = [results_array transpose(mu) transpose(var)]
end

results_array = results_array[:,2:end]

scaled_array = mean_center_array(results_array)
