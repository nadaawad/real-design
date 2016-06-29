`timescale 1ns / 1ps

module top_module(clk,reset,finish,reset_vXv1,reset_mXv1,halt);
	

	parameter element_width = 32;	
	parameter no_of_units = 8;	
	parameter no_of_elements_on_col_nos = 20 ;	 
	parameter no_of_row_by_vector_modules = 4;	  
	parameter no_of_elements_in_p_emap_output = 8;
	parameter col_nos_values_width = 24;
	parameter memory_height = 1000;
	parameter memory_A_height = memory_height * (no_of_units/no_of_row_by_vector_modules) + no_of_units ;
	parameter address_width= $clog2(memory_height)+1;
	parameter A_address_width= $clog2(memory_A_height)+1; 
	
	parameter multiples_memory_value_width=4;
	
	
	
	input wire clk; 
	input wire reset; 
	input reset_vXv1;
	input reset_mXv1;
	output wire finish;
	output wire halt; 
	
	
	//wire write_enable;

	wire [address_width - 1 : 0] memoryX_write_address;
	wire [address_width - 1 : 0] memoryP_write_address;
	wire [address_width - 1 : 0] memoryR_write_address;
	wire [address_width - 1 : 0] memoryX_read_address;	
	
	wire [A_address_width - 1 : 0] memoryA_read_address;	
	wire [A_address_width - 1 : 0] col_nos_read_address;	
	wire [A_address_width - 1 : 0] multiples_read_address;	
	
	wire [address_width - 1 : 0] memoryR_read_address;
	wire [address_width - 1 : 0] memoryP_read_address;
	wire [address_width - 1 : 0] memoryP_v2_read_address;
	
	


	wire [element_width * no_of_units - 1 : 0] memoryP_v2_output;
	wire [no_of_units * element_width - 1 : 0] memoryR_output; 
	
	wire[no_of_row_by_vector_modules*no_of_elements_in_p_emap_output*element_width-1:0] memA_output;
	wire[no_of_row_by_vector_modules*no_of_elements_on_col_nos*col_nos_values_width-1:0] col_nos_output;
	wire [multiples_memory_value_width*no_of_row_by_vector_modules-1:0] multiples_output ;	
	wire [no_of_row_by_vector_modules*no_of_elements_in_p_emap_output*element_width-1:0] Emap_mem_output_row ;
	
	wire [no_of_units * element_width - 1 : 0] memoryX_output;
	wire [no_of_units * element_width - 1 : 0] memoryX_input;
	
	
	wire [no_of_units * element_width - 1 : 0] memoryP_input;
	wire [no_of_units * element_width - 1 : 0] memoryR_input;
	
	wire [no_of_row_by_vector_modules-1:0]you_can_read;
	 wire [no_of_row_by_vector_modules-1:0] I_am_ready;

        
	wire finishP;
	wire finishR;
	
	
	wire result_mem_we_4;

	wire memoryRprev_we;
	
wire result_mem_we_5;
	
wire read_again;

	wire result_mem_we_6;
	
	wire read_again_2;
	
wire[32-1:0] result_mem_counter_5;
	

wire start;
	
wire vXv1_finish;

	wire finish_all;
   
   wire P_Emap_read_preprocess,P_Emap_write_enable ;
   wire[31:0]total_with_additional_A;
   wire[31:0] total;
   wire memories_pre_preprocess;

	
   reg memories_preprocess=0;
   
   
   
   
   main_alu  #(.element_width (element_width),.memory_height(memory_height),.multiples_memory_value_width(multiples_memory_value_width),.address_width(address_width),.no_of_row_by_vector_modules(no_of_row_by_vector_modules))
	alu(total,clk,reset,reset_vXv1,reset_mXv1,memA_output,Emap_mem_output_row,multiples_output,col_nos_output,total_with_additional_A,you_can_read,memories_pre_preprocess,memoryP_v2_output,memoryR_output,memoryR_read_address,memoryX_output,memoryP_input,memoryR_input,memoryX_input,finish,outsider_read_now,result_mem_we_4,memoryRprev_we,result_mem_we_5,result_mem_counter_5,read_again,start,read_again_2,result_mem_we_6,vXv1_finish,finish_all, I_am_ready);
	
	
	
	
	
 control_unit #(.no_of_units(no_of_units),.element_width (element_width ),.address_width(address_width),.A_address_width(A_address_width)) 
	CU(total,clk,reset,finish,memories_pre_preprocess,memoryP_write_enable,memoryR_write_enable,memoryX_write_enable,memoryA_read_address,memoryP_read_address,memoryP_v2_read_address,memoryR_read_address,memoryX_read_address,memoryP_write_address,memoryR_write_address ,memoryX_write_address,halt,reset_vXv1,outsider_read_now,result_mem_we_4,memoryRprev_we,result_mem_we_5,result_mem_counter_5,read_again,start,read_again_2,result_mem_we_6,vXv1_finish,finish_all);
	
	
	
	
	genvar j;
		generate
		for(j=0;j<no_of_row_by_vector_modules;j=j+1) begin:instantiate_P_Emap_8	
			P_Emap_8 #(.no_of_units(no_of_units),.element_width(element_width),.no_of_elements_on_col_nos(no_of_elements_on_col_nos),.no_of_elements_in_output(no_of_elements_in_p_emap_output),
			.memory_height(memory_height),.address_width(address_width),
			.col_nos_values_width(col_nos_values_width),.multiples_memory_value_width(multiples_memory_value_width)) 
			P_Emap_mem (clk, memoryP_input,memories_preprocess,memoryP_write_enable, memoryP_write_address,
			col_nos_output[(no_of_row_by_vector_modules-j)*no_of_elements_on_col_nos*col_nos_values_width-1-:no_of_elements_on_col_nos*col_nos_values_width],
			Emap_mem_output_row[((no_of_row_by_vector_modules-j))*no_of_elements_in_p_emap_output*element_width-1-:(no_of_elements_in_p_emap_output*element_width)],
			multiples_output[(no_of_row_by_vector_modules-j)*multiples_memory_value_width-1-:multiples_memory_value_width],you_can_read[no_of_row_by_vector_modules-j-1],I_am_ready[no_of_row_by_vector_modules-j-1]);
		end
	endgenerate


	memP_v2 #(.element_width (element_width ),.memory_height(memory_height),.address_width(address_width)) 
	matP_v2(clk, memoryP_input, memoryP_write_enable,memoryP_v2_read_address,memoryP_write_address, memoryP_v2_output,finishP_v2);
	
	
	memR #(.element_width (element_width ),.memory_height(memory_height),.address_width(address_width))
	matR(clk,memoryR_input, memoryR_write_enable,memoryR_read_address,memoryR_write_address, memoryR_output,finishR);
	
	memA #(.no_of_elements_on_col_nos(no_of_elements_on_col_nos),.no_of_row_by_vector_modules(no_of_row_by_vector_modules),
	.element_width (element_width ),.memory_A_height(memory_A_height),.address_width(A_address_width),.multiples_memory_value_width(multiples_memory_value_width),
	.no_of_elements_in_output(no_of_elements_in_p_emap_output))
	matA(clk,memoryA_read_address,memA_output,memories_preprocess,multiples_output, I_am_ready);	 
	
	
	col_nos #(.no_of_elements_on_col_nos(no_of_elements_on_col_nos),.no_of_row_by_vector_modules(no_of_row_by_vector_modules),
	.memory_A_height(memory_A_height),.address_width(A_address_width),.col_nos_values_width(col_nos_values_width))
	col_nos_memory(clk,memoryA_read_address,col_nos_output);
	
	multiples_memory #(.no_of_row_by_vector_modules(no_of_row_by_vector_modules),
	.memory_A_height(memory_A_height),.address_width(A_address_width),.multiples_memory_value_width(multiples_memory_value_width))
	multiples_mat (clk,memoryA_read_address,multiples_output);
	
	memX #(.element_width (element_width ),.memory_height(memory_height),.address_width(address_width))
	matX(clk, memoryX_input, memoryX_write_enable,memoryX_read_address,memoryX_write_address, memoryX_output);
	
 	parameters_mem parameters_memory(total_with_additional_A,total);

      

	always @(posedge clk)
		begin
		 memories_preprocess <= memories_pre_preprocess	;	
		   
			
		end	



        
        
   	
	
endmodule
