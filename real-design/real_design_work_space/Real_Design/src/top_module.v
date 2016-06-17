`timescale 1ns / 1ps

module top_module(clk,reset,finish,reset_vXv1,reset_mXv1,halt);
	
	parameter number_of_clusters =1;
	parameter number_of_equations_per_cluster =19;
	parameter element_width = 32;
	parameter memories_address_width=20;// m7tag t3deel	
	parameter no_of_units = 8;
	
	input wire clk; 
	input wire reset; 
	input reset_vXv1;
	input reset_mXv1;
	
	output wire halt; 
	output wire finish;
	
	
	wire [memories_address_width - 1 : 0] memoryX_write_address;
	wire [memories_address_width - 1 : 0] memoryP_write_address; 
	wire [memories_address_width - 1 : 0] memoryP_write_address_old;
	wire [memories_address_width - 1 : 0] memoryR_write_address;
	wire [memories_address_width - 1 : 0] memoryX_read_address;
	wire [memories_address_width - 1 : 0] memoryA_read_address;
	wire [memories_address_width - 1 : 0] memoryR_read_address;
	wire [memories_address_width - 1 : 0] memoryP_read_address;
	wire [memories_address_width - 1 : 0] memoryP_v2_read_address;
	
	

	wire [element_width * number_of_equations_per_cluster - 1 : 0] memoryP_output;
	wire [element_width * no_of_units - 1 : 0] memoryP_v2_output;
	wire [element_width * no_of_units - 1 : 0] memoryR_output;
	wire [element_width*(3* number_of_equations_per_cluster-2*2+2)-1 : 0] memoryA_output;
	wire [number_of_equations_per_cluster * element_width - 1 : 0] memoryX_output;
	wire [no_of_units * element_width - 1 : 0] memoryX_input;
	wire [no_of_units * element_width - 1 : 0] memoryP_input;
	wire [no_of_units * element_width - 1 : 0] memoryP_input_old;
	wire [number_of_equations_per_cluster * element_width - 1 : 0] memoryR_input;
	


        
	wire finishP;
	wire finishR;
	
	wire result_mem_we_4;
    wire memoryRprev_we;   
	wire result_mem_we_5;
	wire read_again;
	wire result_mem_we_6;
	wire read_again_2;
	wire[31:0] result_mem_counter_5;
    
	wire start;	
	wire vXv1_finish;
	wire finish_all;


	main_alu #(.number_of_clusters(number_of_clusters),.number_of_equations_per_cluster(number_of_equations_per_cluster),.element_width (element_width ))
	alu(clk,reset,reset_vXv1,reset_mXv1,memoryA_output,memoryP_output,memoryP_v2_output,memoryR_output,memoryR_read_address,memoryX_output,memoryP_input,memoryR_input,memoryX_input,finish,mXv1_finish,result_mem_we_4,memoryRprev_we,result_mem_we_5,result_mem_counter_5,read_again,start,read_again_2,result_mem_we_6,vXv1_finish,finish_all);
	
	control_unit #(.no_of_units(no_of_units),.number_of_clusters(number_of_clusters),.number_of_equations_per_cluster(number_of_equations_per_cluster),.element_width (element_width ),.memories_address_width(memories_address_width)) 
	CU(clk,reset,finish,memoryP_write_enable,memoryR_write_enable,memoryX_write_enable,memoryA_read_address,memoryP_read_address,memoryP_v2_read_address,memoryR_read_address,memoryX_read_address,memoryP_write_address,memoryR_write_address ,memoryX_write_address,halt,reset_vXv1,mXv1_finish,result_mem_we_4,memoryRprev_we,result_mem_we_5,result_mem_counter_5,read_again,start,read_again_2,result_mem_we_6,vXv1_finish,finish_all);
	
	memP #(.number_of_clusters(number_of_clusters),.number_of_equations_per_cluster(number_of_equations_per_cluster),.element_width (element_width ),.memories_address_width(memories_address_width)) 
	matP(clk, memoryP_input_old, memoryP_write_enable_old,memoryP_read_address,memoryP_write_address_old, memoryP_output,finishP);
	
	memP_v2 #(.number_of_clusters(number_of_clusters),.number_of_equations_per_cluster(number_of_equations_per_cluster),.element_width (element_width ),.memories_address_width(memories_address_width)) 
	matP_v2(clk, memoryP_input, memoryP_write_enable,memoryP_v2_read_address,memoryP_write_address, memoryP_v2_output,finishP_v2);
	
	
	memR #(.number_of_clusters(number_of_clusters),.number_of_equations_per_cluster(number_of_equations_per_cluster),.element_width (element_width ))
	matR(clk,memoryR_input, memoryR_write_enable,memoryR_read_address,memoryR_write_address, memoryR_output,finishR);
	
	memA #(.number_of_clusters(number_of_clusters),.number_of_equations_per_cluster(number_of_equations_per_cluster),.element_width (element_width ),.memories_address_width(memories_address_width))
	matA(memoryA_read_address,clk,memoryA_output);
	
	memX #(.number_of_clusters(number_of_clusters),.number_of_equations_per_cluster(number_of_equations_per_cluster),.element_width (element_width ),.memories_address_width(memories_address_width))
	matX(clk, memoryX_input, memoryX_write_enable,memoryX_read_address,memoryX_write_address, memoryX_output);


      




        
        
   	
	
endmodule
