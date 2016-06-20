`timescale 1ns / 1ps

module top_module(clk,reset,reset_vXv1,reset_mXv1,halt,reset_cluster);
	

	parameter element_width = 32;
	parameter memories_address_width=32;	
	parameter no_of_units = 8;	
	parameter no_of_elements_on_col_nos = 20 ;	 
	parameter no_of_row_by_vector_modules = 4;	  
	parameter no_of_elements_in_p_emap_output = 8;
	
	input wire clk; 
	input wire reset; 
	input reset_vXv1;
	input reset_mXv1;
	
	output wire halt; 
	output wire reset_cluster;
	
	//wire write_enable;
	wire finish;
	wire [memories_address_width - 1 : 0] memoryX_write_address;
	wire [memories_address_width - 1 : 0] memoryP_write_address;
	wire [memories_address_width - 1 : 0] memoryR_write_address;
	wire [memories_address_width - 1 : 0] memoryX_read_address;	
	
	wire [memories_address_width - 1 : 0] memoryA_read_address;	
	wire [memories_address_width - 1 : 0] col_nos_read_address;	
	wire [memories_address_width - 1 : 0] multiples_read_address;	
	
	wire [memories_address_width - 1 : 0] memoryR_read_address;
	wire [memories_address_width - 1 : 0] memoryP_read_address;
	wire [memories_address_width - 1 : 0] memoryP_v2_read_address;
	
	


	wire [element_width * no_of_units - 1 : 0] memoryP_v2_output;
	wire [no_of_units * element_width - 1 : 0] memoryR_output; 
	
	wire[no_of_row_by_vector_modules*no_of_elements_on_col_nos*element_width-1:0] memA_output;
	wire[no_of_row_by_vector_modules*no_of_elements_on_col_nos*32-1:0] col_nos_output;
	wire [32*no_of_row_by_vector_modules-1:0] multiples_output ;	
	wire [no_of_row_by_vector_modules*no_of_elements_in_p_emap_output*element_width-1:0] Emap_mem_output_row ;
	
	wire [no_of_units * element_width - 1 : 0] memoryX_output;
	wire [no_of_units * element_width - 1 : 0] memoryX_input;
	wire [no_of_units * element_width - 1 : 0] memoryR_input;
	


        //wire finish_memory;
	wire finishP;
	wire finishR;
	//wire finishA;
	//wire finishX;
	wire finish_alpha;
	wire finish_iteration;
	wire flag_in_1;
	wire flag_in_2;
	wire flag_out_1;
	wire flag_out_2;
	wire flag2;	  
	wire result_mem_we_4;
    wire memoryRprev_we;   
	wire result_mem_we_5;
	wire read_again;
	wire[31:0] result_mem_counter_5;
   wire[element_width*no_of_units-1:0] result_mem_in_5;
   
   wire P_Emap_read_preprocess,P_Emap_write_enable ;
   wire[31:0]total_with_additional_A;


	main_alu #(.element_width (element_width ))
	alu(clk,reset,reset_vXv1,reset_mXv1,finish_alpha,memoryA_output,memoryP_output,memoryP_v2_output,memoryR_output,memoryR_read_address,memoryX_output,memoryP_input,memoryR_input,memoryX_input,finish,finish_iteration,outsider_read_now,result_mem_we_4,memoryRprev_we,result_mem_we_5,result_mem_counter_5,result_mem_in_5,read_again);
	
	control_unit #(.no_of_units(no_of_units),.element_width (element_width ),.memories_address_width(memories_address_width)) 
	CU(clk,reset,finish,finish_alpha,memoryP_write_enable,memoryR_write_enable,memoryX_write_enable,memoryA_read_address,memoryP_read_address,memoryP_v2_read_address,memoryR_read_address,memoryX_read_address,memoryP_write_address,memoryR_write_address ,memoryX_write_address,halt,reset_cluster,finish_iteration,reset_vXv1,outsider_read_now,result_mem_we_4,memoryRprev_we,result_mem_we_5,result_mem_counter_5,read_again);

	
		genvar j;
	generate
	for(j=0;j<no_of_row_by_vector_modules;j=j+1) begin:instantiate_P_Emap_8
		
	P_Emap_8 #(.no_of_units(no_of_units),.element_width(element_width),.no_of_elements_on_col_nos(no_of_elements_on_col_nos),.no_of_elements_in_output(no_of_elements_in_p_emap_output)) 
	P_Emap_mem (clk,P_Emap_read_preprocess,P_Emap_write_enable,
	col_nos_output[(no_of_row_by_vector_modules-j)*no_of_elements_on_col_nos*32-1-:no_of_elements_on_col_nos*32],
	Emap_mem_output_row[((no_of_row_by_vector_modules-j))*no_of_elements_in_p_emap_output*element_width-1-:(no_of_elements_in_p_emap_output*element_width)],
	multiples_output[(no_of_row_by_vector_modules-j)*32-1-:32]);
	end
	endgenerate


	memP_v2 #(.element_width (element_width ),.memories_address_width(memories_address_width)) 
	matP_v2(clk, memoryP_input, memoryP_write_enable,memoryP_v2_read_address,memoryP_write_address, memoryP_v2_output,finishP_v2);
	
	
	memR #(.element_width (element_width ))
	matR(clk,result_mem_in_5, memoryR_write_enable,memoryR_read_address,memoryR_write_address, memoryR_output,finishR);
	
	memA #(.no_of_elements_on_col_nos(no_of_elements_on_col_nos),.no_of_row_by_vector_modules(no_of_row_by_vector_modules),.element_width (element_width ))
	matA(clk,memoryA_read_address,memA_output);	 
	
	col_nos #(.no_of_elements_on_col_nos(no_of_elements_on_col_nos),.no_of_row_by_vector_modules(no_of_row_by_vector_modules))
	col_nos_memory(clk,col_nos_read_address,col_nos_output);
	
	multiples_memory #(.no_of_row_by_vector_modules(no_of_row_by_vector_modules))
	multiples_mat (clk,multiples_read_address,multiples_output);
	
	memX #(.element_width (element_width ),.memories_address_width(memories_address_width))
	matX(clk, memoryX_input, memoryX_write_enable,memoryX_read_address,memoryX_write_address, memoryX_output);
	
 	parameters_mem parameters_memory(total_with_additional_A);

      




        
        
   	
	
endmodule
