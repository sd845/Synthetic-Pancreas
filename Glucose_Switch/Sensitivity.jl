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
    time_stop = 16.0
    time_step_size = 0.01
    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"
    #path_to_data_dir = "$(pwd())/data"
    # Load the data dictionary (uses the default biophysical_constants file)
    model_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)
    # Phase 1: parameter update =========================================================================== #
    # update the paramaters in the model data dictionary -
    # for now - lets only search over dG's -
    R = model_data_dictionary["R"]
    T_K = model_data_dictionary["T_K"]
    # compute W -
    tmp_W_array = Float64[]
    for index = 1:13
        parameter_guess = parameter_guess_array[index]
        value = exp(-1*parameter_guess/(R*T_K))
        push!(tmp_W_array,value)
    end
    # update the control W's -
    control_parameter_dictionary = model_data_dictionary["control_parameter_dictionary"]
    control_parameter_dictionary["W_GFPssrA_RNAP"] = tmp_W_array[1]
    control_parameter_dictionary["W_GFPssrA_sigma_70"] = tmp_W_array[2]
    control_parameter_dictionary["W_GFPssrA_GntR"] = tmp_W_array[3]
    control_parameter_dictionary["W_GntR_RNAP"] = tmp_W_array[4]
    control_parameter_dictionary["W_GntR_sigma_70"] = tmp_W_array[5]
    control_parameter_dictionary["W_RFPssrA_RNAP"] = tmp_W_array[6]
    control_parameter_dictionary["W_RFPssrA_sigma_28"] = tmp_W_array[7]
    control_parameter_dictionary["W_RFPssrA_aplha_sigma28"] = tmp_W_array[8]
    control_parameter_dictionary["W_aplha_sigma28_RNAP"] = tmp_W_array[9]
    control_parameter_dictionary["W_aplha_sigma28_sigma_70"] = tmp_W_array[10]
    control_parameter_dictionary["W_aplha_sigma28_GntR"] = tmp_W_array[11]
    control_parameter_dictionary["W_sigma_28_RNAP"] = tmp_W_array[12]
    control_parameter_dictionary["W_sigma_28_sigma_70"] = tmp_W_array[13]
    model_data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
    binding_parameter_dictionary = model_data_dictionary["binding_parameter_dictionary"]
    binding_parameter_dictionary["n_GFPssrA_sigma_70"] = parameter_guess_array[14]
    binding_parameter_dictionary["K_GFPssrA_sigma_70"] = parameter_guess_array[15]
    binding_parameter_dictionary["n_GFPssrA_GntR"] = parameter_guess_array[16]
    binding_parameter_dictionary["K_GFPssrA_GntR"] = parameter_guess_array[17]
    binding_parameter_dictionary["n_GntR_sigma_70"] = parameter_guess_array[18]
    binding_parameter_dictionary["K_GntR_sigma_70"] = parameter_guess_array[19]
    binding_parameter_dictionary["n_RFPssrA_sigma_28"] = parameter_guess_array[20]
    binding_parameter_dictionary["K_RFPssrA_sigma_28"] = parameter_guess_array[21]
    binding_parameter_dictionary["n_RFPssrA_aplha_sigma28"] = parameter_guess_array[22]
    binding_parameter_dictionary["K_RFPssrA_aplha_sigma28"] = parameter_guess_array[23]
    binding_parameter_dictionary["n_aplha_sigma28_sigma_70"] = parameter_guess_array[24]
    binding_parameter_dictionary["K_aplha_sigma28_sigma_70"] = parameter_guess_array[25]
    binding_parameter_dictionary["n_aplha_sigma28_GntR"] = parameter_guess_array[26]
    binding_parameter_dictionary["K_aplha_sigma28_GntR"] = parameter_guess_array[27]
    binding_parameter_dictionary["n_sigma_28_sigma_70"] = parameter_guess_array[28]
    binding_parameter_dictionary["K_sigma_28_sigma_70"] = parameter_guess_array[29]
    model_data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
    # time constant modifier -
    time_constant_modifier_array = [
        0.0 ;   # 1 GFPssrA
        0.0 ;   # 2 GntR
        0.0 ;   # 3 RFPssrA
        0.0 ;   # 4 aplha_sigma28
        0.0 ;   # 5 sigma_28
        0.0 ;   # 6 sigma_70
        parameter_guess_array[30]   ;   # 7 mRNA_GFPssrA
        parameter_guess_array[31]   ;   # 8 mRNA_GntR
        parameter_guess_array[32]   ;   # 9 mRNA_RFPssrA
        parameter_guess_array[33]   ;   # 10    mRNA_aplha_sigma28
        parameter_guess_array[34]   ;   # 11    mRNA_sigma_28
        1.0 ;                           # 12    mRNA_sigma_70
        parameter_guess_array[35]   ;   # 13    protein_GFPssrA
        parameter_guess_array[36]   ;   # 14    protein_GntR
        parameter_guess_array[37]   ;   # 15    protein_RFPssrA
        parameter_guess_array[38]   ;   # 16    protein_aplha_sigma28
        parameter_guess_array[39]   ;   # 17    protein_sigma_28
        1.0 ;                           # 18    protein_sigma_70
        ]
    model_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array
    # setup degradation_modifier_array -
    degradation_modifier_array = [
        0.0                         ;   # 1     GFPssrA
        0.0                         ;   # 2     GntR
        0.0                         ;   # 3     RFPssrA
        0.0                         ;   # 4     aplha_sigma28
        0.0                         ;   # 5     sigma_28
        0.0                         ;   # 6     sigma_70
        parameter_guess_array[40]   ;   # 7     mRNA_GFPssrA
        parameter_guess_array[41]   ;   # 8     mRNA_GntR
        parameter_guess_array[42]   ;   # 9     mRNA_RFPssrA
        parameter_guess_array[43]   ;   # 10    mRNA_aplha_sigma28
        parameter_guess_array[44]   ;   # 11    mRNA_sigma_28
        1.0                         ;   # 12    mRNA_sigma_70
        parameter_guess_array[45]   ;   # 13    protein_GFPssrA
        parameter_guess_array[46]   ;   # 14    protein_GntR
        parameter_guess_array[47]   ;   # 15    protein_RFPssrA
        parameter_guess_array[48]   ;   # 16    protein_aplha_sigma28
        parameter_guess_array[49]   ;   # 17    protein_sigma_28
        parameter_guess_array[50]   ;   # 18    protein_sigma_70
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
    # ===================================================================================================== #
    # Phase 3: compute the model performance metrics ====================================================== #
    p_GFP_AUC = integrate(TSIM,XSIM[:,index])
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
            0 1        ;   # 1 W_GFPssrA_RNAP
            0 1   ;   # 2 W_deGFPssrA_sigma_70
            0 1    ;   # 3 W_GFPssrA_GntR
            0 1       ;   # 4 W_GntR_RNAP
            0 1   ;   # 5 W_GntR_sigma_70
            0 1        ;   # 6 W_RFPssrA_RNAP
            0 1   ;   # 7 W_RFPssrA_sigma_28
            0 1   ;   # 8 W_RFPssrA_aplha_sigma28
            0 1     ;   # 9 W_aplha_sigma28_RNAP
            0 1  ;   # 10 W_aplha_sigma28_sigma_70
            0 1   ;   # 11 W_aplha_sigma28_GntR
            0 1    ;   # 12 W_sigma_28_RNAP
            0 1   ;   # 13 W_sigma_28_sigma_70
            # W_sigma_70_RNAP = 0.0
            # binding parameters -
            0.50 2.2   		; #14  n_GFPssrA_sigma_70
            26.0 100.0 		; #15  K_GFPssrA_sigma_70
            1.5 3.2 		; #16  n_GFPssrA_GntR
            0.75 1.25		; #17  K_GFPssrA_GntR
            1.3 1.8			; #18  n_GntR_sigma_70
            62.0 100.0  	; #19  K_GntR_sigma_70
            1.3 1.8			; #20  n_RFPssrA_sigma_28
            62.0 100.0  	; #21  K_RFPssrA_sigma_28
            1.3 1.8			; #22  n_RFPssrA_aplha_sigma28
            62.0 100.0  	; #23  K_RFPssrA_aplha_sigma28
            1.3 1.8			; #24  n_aplha_sigma28_sigma_70
            62.0 100.0  	; #25  K_aplha_sigma28_sigma_70
            1.3 1.8			; #26  n_aplha_sigma28_GntR
            62.0 100.0  	; #27  K_aplha_sigma28_GntR
            1.3 1.8			; #28  n_sigma_28_sigma_70
            62.0 100.0  	; #29  K_sigma_28_sigma_70
            # time constants -
            0.05 0.067      ;   # 30    mRNA_GFPssrA
            0.05 0.067      ;   # 31    mRNA_GntR
            0.05 0.067      ;   # 32    mRNA_RFPssrA
            0.05 0.067      ;   # 33    mRNA_aplha_sigma28
            0.05 0.067      ;   # 34    mRNA_sigma_28
            0.45 0.7        ;   # 35    protein_GFPssrA
            0.45 0.7        ;   # 36    protein_GntR
            0.45 0.7        ;   # 37    protein_RFPssrA
            0.45 0.7        ;   # 38    protein_aplha_sigma28
            0.45 0.7        ;   # 39    protein_sigma_28
            # mRNA_sigma_70 = 1.0 and protein_sigma_70 = 1.0
            # degradation mods -
            1.67 3.07	;	# 40	mRNA_GFPssrA
            1.67 3.07	;	# 41	mRNA_GntR
            1.67 3.07	;	# 42	mRNA_RFPssrA
            1.67 3.07	;	# 43	mRNA_aplha_sigma28
            1.67 3.07	;	# 44	mRNA_sigma_28
            50 60 		;	# 45	protein_GFPssrA
            50 60 		;	# 46	protein_GntR
            50 60 		;	# 47	protein_RFPssrA
            50 60 		;	# 48	protein_aplha_sigma28
            50 60 		;	# 49	protein_sigma_28
            0.15 0.32 	;	# 50	protein_sigma_70
        ];
    #@show index
    number_of_parameters = size(sample_bounds_array)[1]
            for parameter_index = 1:(number_of_parameters)
                # get row of parameters -
                lb = sample_bounds_array[parameter_index,1]
                ub = sample_bounds_array[parameter_index,2]
                # create the tuple -
                tmp_tuple = (lb,ub)
                #display(tmp_tuple)
                # cache -
                push!(parameter_bounds_array,tmp_tuple)
            end
    # do the global sensitivity analysis -
    sensitivity_results = gsa(SF,Morris(total_num_trajectory=10000,num_trajectory=1000),parameter_bounds_array)
    # return -
    return sensitivity_results
end
# setup paths -
#path_to_ensemble_file = "$(pwd())/Ensemble-c1-restriction-T20.dat"
# compute a sensitivity array for the AUC of each species -
species_index_array = [7 8 9 10 11 13 14 15 16 17]
#species_index_array = [7]
number_of_species = length(species_index_array)
number_of_parameters = 50
results_array = zeros(number_of_parameters,1)
for species_index in species_index_array
    @show species_index
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
#scaled_array = mean_center_array(results_array)
