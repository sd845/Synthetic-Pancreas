using DelimitedFiles
using Statistics
using PyPlot
using PyCall
@pyimport matplotlib.patches as patches

# define my colors -
negative_color = (1/255)*[255,51,51]
positive_color = (1/255)*[51,153,255]

function calculate_x_axis_tick_positions(data_array;box_width=1.0)

    # what is the size?
    (number_of_rows, number_of_cols) = size(data_array)

    # initialize -
    x_axis_points = Float64[]
    epsilon = 0.2;

    for col_index = 1:number_of_cols

        # compute origin point -
        origin_point = (col_index - 1)+(col_index - 1)*epsilon + 1
        push!(x_axis_points,origin_point+0.5*box_width)
    end

    return x_axis_points
end

function calculate_y_axis_tick_positions(data_array; box_height=1.0)

    # what is the size?
    (number_of_rows, number_of_cols) = size(data_array)

    # initialize -
    y_axis_points = Float64[]

    for col_index = 1:number_of_cols

        # get the data array -
        data_scaled = data_array[:,col_index]

        # how many patches per col?
        number_of_patches = length(data_scaled)
        epsilon = 0.2;
        for row_index = 1:number_of_patches

            # compute origin point -
            origin_point = [(col_index - 1)+(col_index - 1)*epsilon + 1,(row_index - 1)+(row_index - 1)*epsilon+1];

            push!(y_axis_points,origin_point[2]+0.5*box_height)
        end
    end

    return y_axis_points
end

function visualize(data_array,x_axis_ticks,y_axis_ticks,y_label_array)

    # what is the size?
    (number_of_rows, number_of_cols) = size(data_array)

    epsilon = 0.2;
    ax = gca()

    # set the axis tick postions -
    #ax.set_xticks(x_axis_points)
    ax.set_yticks(y_axis_ticks)
    ax.set_yticklabels(y_label_array, fontsize=7,fontname="Arial")
    subplots_adjust(left=0.25)

    # create x-label array -
    x_label_array = ["S$(x)" for x=1:number_of_cols]
    #ax.set_xticklabels(x_label_array, rotation=45,fontsize=8)
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.spines["bottom"].set_visible(false)
    #ax.spines["left"].set_visible(false)

    ax.get_xaxis().set_visible(false)

    for col_index = 1:number_of_cols

        # get the data array -
        data_scaled = data_array[:,col_index]

        # how many patches per col?
        number_of_patches = length(data_scaled)
        for row_index = 1:number_of_patches

            # compute origin point -
            origin_point = [(col_index - 1)+(col_index - 1)*epsilon + 1,(row_index - 1)+(row_index - 1)*epsilon+1];

            # compute alpha -
            data_value = data_array[row_index, col_index]
            if (data_value == 0.0)

                # compute a new color -
                patch_face_color = "#ffffff"

                # draw the square -
                ax.add_patch(

                           patches.Rectangle(origin_point,   # (x,y)
                               1.0,          # width
                               1.0,          # height
                               facecolor=patch_face_color,
                               edgecolor="grey",
                               linewidth=0.5,
                           )
                       )

            else

                #beta = (data_value - a)/(b - a)

                # compute a new color -
                #patch_face_color = (1-beta)*negative_color .+ beta*positive_color
                patch_face_color = "#ffffff"
                if (data_value == 1.0)
                    patch_face_color="#DDE3ED"
                elseif (data_value == 2.0)
                    patch_face_color="#C8D1E0"
                elseif (data_value == 3.0)
                    patch_face_color="#AFBACC"
                elseif (data_value == 4.0)
                    patch_face_color="#8E99AB"
                elseif (data_value == 5.0)
                    patch_face_color="#707A8A"
                elseif (data_value == 6.0)
                    patch_face_color="#58606E"
                elseif (data_value == 7.0)
                    patch_face_color="#434A54"
                elseif (data_value == 8.0)
                    patch_face_color="#333840"
                end

                # draw the square -
                ax.add_patch(

                           patches.Rectangle(origin_point,   # (x,y)
                               1.0,          # width
                               1.0,          # height
                               facecolor=patch_face_color,
                               edgecolor="grey",
                               linewidth=0.5,
                               alpha = 1.0
                           )
                       )
            end
        end
    end


    axis("square")
end

function reduce(path_to_sensitivity_array::String; epsilon::Float64=1e-5)

    # parameter name array =
    parameter_name_array = [
		"dG_GFPssrA_RNAP"	            ;	# 1
		"dG_deGFPssrA_sigma_70"			; 	# 2
		"dG_GFPssrA_GntR"				;	# 3
		"dG_GntR_RNAP"					;	# 4
		"dG_GntR_sigma_70"	            ;	# 5
		"dG_RFPssrA_RNAP"	            ;	# 6
		"dG_RFPssrA_sigma_28"	        ;	# 7
		"dG_RFPssrA_aplha_sigma28"	    ;	# 8
		"dG_aplha_sigma28_RNAP"	        ;	# 9
		"dG_aplha_sigma28_sigma_70"	    ;	# 10
		"dG_aplha_sigma28_GntR"	        ;	# 11
		"dG_sigma_28_RNAP"	            ;	# 12
		"dG_sigma_28_sigma_70"	        ;	# 13
 		"n_GFPssrA_sigma_70"			;	# 14
   		"K_GFPssrA_sigma_70" 			; 	# 15
		"n_GFPssrA_GntR"         		; 	# 16
		"K_GFPssrA_GntR"        		; 	# 17
		"n_GntR_sigma_70"       		;  	# 18
		 "K_GntR_sigma_70"       		; 	# 19
		"n_RFPssrA_sigma_28"          	; 	# 20
		"K_RFPssrA_sigma_28"       		; 	# 21
		"n_RFPssrA_aplha_sigma28"       ; 	# 22
	 	"K_RFPssrA_aplha_sigma28"      	; 	# 23
		"n_aplha_sigma28_sigma_70"      ; 	# 24
		"K_aplha_sigma28_sigma_70"      ; 	# 25
		"n_aplha_sigma28_GntR"        	; 	# 26
		"K_aplha_sigma28_GntR"       	; 	# 27
		"n_sigma_28_sigma_70"        	; 	# 28
		"K_sigma_28_sigma_70"      		; 	# 29
		"tc_mRNA_GFPssrA" 				;	# 30
		"tc_mRNA_GntR"					;	# 31
		"tc_mRNA_RFPssrA"				;	# 32
		"tc_mRNA_aplha_sigma28"			;	# 33
		"tc_mRNA_sigma_28" 				;	# 34
		"tc_protein_GFPssrA"  			;	# 35
		"tc_protein_GntR"				;	# 36
		"tc_protein_RFPssrA" 			;	# 37
		"tc_protein_aplha_sigma28" 		;	# 38
		"tc_protein_sigma_28" 			;	# 39
		"deg_mRNA_GFPssrA"				;	# 40
		"deg_mRNA_GntR"					;	# 41
		"deg_mRNA_RFPssrA"				;	# 42
		"deg_mRNA_aplha_sigma28"		;	# 43
		"deg_mRNA_sigma_28"				;	# 44
		"deg_protein_GFPssrA"			;	# 45
		"deg_protein_GntR"				;	# 46
		"deg_protein_RFPssrA"			;	# 47
		"deg_protein_aplha_sigma28"		;	# 48
		"deg_protein_sigma_28"			;	# 49
		"deg_protein_sigma_70"			;	# 50
    ];

    # load -
    sensitivity_array = readdlm(path_to_sensitivity_array)

    # parameters to kepp -
    index_keep_p = Int64[]

    # the parameters are on the rows - are any values abpove a threshold -
    (NR,NC) = size(sensitivity_array)
    for row_index = 1:NR

        # get the row data -
        row_data = sensitivity_array[row_index, :]

        # check -
        idx_check = findall(x->x>=epsilon, row_data)
        if (isempty(idx_check) == false)
            push!(index_keep_p, row_index)
        end
    end

    # pull out rows -
    reduced_array = sensitivity_array[index_keep_p, :]

    # last thing - lets put things in groups of 0,1,2,3 (or none,low,medium,high)
    (NR,NC) = size(reduced_array)
    category_array = zeros(NR,NC)
    for row_index = 1:NR
        for col_index = 1:NC

            old_value = reduced_array[row_index, col_index]
            order_of_magnitude = floor(Int,log10(old_value))
            oom_eps = floor(Int,log10(epsilon))

            # rules -
            if (order_of_magnitude<-4)                                  # below -4
                category_array[row_index, col_index] = 0
            elseif (order_of_magnitude == -4)                           # O(-3)
                category_array[row_index, col_index] = 1
            elseif (order_of_magnitude == -3)                           # O(-2)
                category_array[row_index, col_index] = 2
            elseif (order_of_magnitude == -2)                           # O(-1)
                category_array[row_index, col_index] = 3
            elseif (order_of_magnitude == -1 )                           # O(0)
                category_array[row_index, col_index] = 4
            elseif (order_of_magnitude == 0)                            # O(1)
                category_array[row_index, col_index] = 5
            elseif (order_of_magnitude == 1)                            # O(2)
                category_array[row_index, col_index] = 6
            elseif (order_of_magnitude == 2)                            # O(3)
                category_array[row_index, col_index] = 7
            else
                category_array[row_index, col_index] = 8
            end
        end
    end

    # grab p labels -
    p_label_array = parameter_name_array[index_keep_p]

    # return -
    return (index_keep_p, p_label_array, category_array)
end

# setup path to sensitvity array -
path_to_sensitivity_array = "./sensitivity_res_mod.dat"

# execute -
(idx_p, p_label_array, rsa) = reduce(path_to_sensitivity_array; epsilon=1e-1)

# calculate the tick locations -
x_axis_points = calculate_x_axis_tick_positions(rsa)
y_axis_points = calculate_y_axis_tick_positions(rsa)

# call plotting code -
#using Plots
#pyplot()
using PlotlyJS
clf()
visualize(rsa,x_axis_points,y_axis_points, p_label_array)
PyPlot.savefig("Sensitivity_results_y.pdf")
#png("plot3")
