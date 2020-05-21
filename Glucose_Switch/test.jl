parameter_bounds_array = Array{Tuple,1}()
sample_bounds_array = [

        # dG's -
        39995.69565886529 68282.2841495398      ;   # 1 W_GFPssrA_RNAP
        -35689.208748865043 -24993.98186199351  ;   # 2 W_deGFPssrA_sigma_70
    -50000.0 25000.0   			  			;   # 3 W_GFPssrA_GntR
    39995.69565886529 68282.2841495398      ;   # 4 W_GntR_RNAP
    -50000.0 25000.0   			  			;   # 5 W_GntR_sigma_70
        39995.69565886529 68282.2841495398      ;   # 6 W_RFPssrA_RNAP
    -50000.0 25000.0   			 		    ;   # 7 W_RFPssrA_sigma_28
    -50000.0 25000.0   			  			;   # 8 W_RFPssrA_aplha_sigma28
        39995.69565886529 68282.2841495398      ;   # 9 W_aplha_sigma28_RNAP
    -50000.0 25000.0   			 			;   # 10 W_aplha_sigma28_sigma_70
    -50000.0 25000.0   			 		    ;   # 11 W_aplha_sigma28_GntR
    39995.69565886529 68282.2841495398      ;   # 12 W_sigma_28_RNAP
    -50000.0 25000.0   			  			;   # 13 W_sigma_28_sigma_70
    # W_sigma_70_RNAP = 0.0


        # binding parameters -
    0.5037004212911795 2.2   	; #14  n_GFPssrA_sigma_70
    25.539812443045424 100.0 	; #15  K_GFPssrA_sigma_70
    0.5 10.0          			; #16  n_GFPssrA_GntR
    0.001 100.0         		; #17  K_GFPssrA_GntR
    0.5 10.0         			; #18  n_GntR_sigma_70
    0.001 100.0         		; #19  K_GntR_sigma_70
    0.5 10.0          			; #20  n_RFPssrA_sigma_28
    0.001 100.0         		; #21  K_RFPssrA_sigma_28
    0.5 10.0         			; #22  n_RFPssrA_aplha_sigma28
    0.001 100.0        			; #23  K_RFPssrA_aplha_sigma28
    0.5 10.0         	 		; #24  n_aplha_sigma28_sigma_70
    0.001 100.0         		; #25  K_aplha_sigma28_sigma_70
    0.5 10.0         	 		; #26  n_aplha_sigma28_GntR
    0.001 100.0        			; #27  K_aplha_sigma28_GntR
    0.5 10.0         			; #28  n_sigma_28_sigma_70
    0.001 100.0        			; #29  K_sigma_28_sigma_70


        # time constants -
    0.0039273045649573235 7.17667775643900 	;	# 30	mRNA_GFPssrA
    0.0001 100.0 	;	# 31	mRNA_GntR
    0.0001 100.0 	;	# 32	mRNA_RFPssrA
    0.0001 100.0 	;	# 33	mRNA_aplha_sigma28
    0.0001 100.0 	;	# 34	mRNA_sigma_28
    0.35035678021027006 10.203484322453043  	;	# 35	protein_GFPssrA
    0.0001 100.0 	;	# 36	protein_GntR
    0.0001 100.0 	;	# 37	protein_RFPssrA
    0.0001 100.0 	;	# 38	protein_aplha_sigma28
    0.0001 100.0 	;	# 39	protein_sigma_28
    # mRNA_sigma_70 = 1.0 and protein_sigma_70 = 1.0

        # degradation mods -
    0.1 10.0	;	# 40	mRNA_GFPssrA
    0.1 10.0	;	# 41	mRNA_GntR
    0.1 10.0	;	# 42	mRNA_RFPssrA
    0.1 10.0	;	# 43	mRNA_aplha_sigma28
    0.1 10.0	;	# 44	mRNA_sigma_28
    0.1 10.0	;	# 45	mRNA_sigma_70
    40.0 60.0	;	# 46	protein_GFPssrA
    0.1 10.0	;	# 47	protein_GntR
    0.1 10.0	;	# 48	protein_RFPssrA
    0.1 10.0	;	# 49	protein_aplha_sigma28
    0.1 10.0	;	# 50	protein_sigma_28
    0.1 10.0	;	# 51	protein_sigma_70

    ];

  number_of_parameters = size(sample_bounds_array)[1]

  for parameter_index = 1:(number_of_parameterssample_bounds_array)
    # get row of parameters -
    lb = sample_bounds_array[parameter_index,1]
    ub = sample_bounds_array[parameter_index,2]

    # create the tuple -
    #tmp_tuple = (lb,ub)

    # cache -
    push!(parameter_bounds_array,(lb,ub))
  end
