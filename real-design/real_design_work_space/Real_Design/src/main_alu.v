
module main_alu(clk,reset,reset_vXv1,reset_mXv1,memoryA_output,memoryP_output,pKold_v2,memoryR_output,memoryR_read_address,memoryX_output,memoryP_input,memoryR_input,memoryX_input,finish,mXv1_finish,result_mem_we_4,memoryRprev_we,result_mem_we_5,result_mem_counter_5,read_again,start,read_again_2,result_mem_we_6,vXv1_finish ,finish_all
	);

    parameter number_of_clusters =40;
	parameter number_of_equations_per_cluster =19;
	parameter element_width = 32;
	parameter memories_address_width=20;   
	parameter no_of_units = 8; 
	parameter memory_read_address_width=20;	

	
	input wire clk;
    input wire reset;
	
    input wire reset_vXv1;
    input wire reset_mXv1; 
    input wire [element_width*(3* number_of_equations_per_cluster-2*2+2)-1 : 0] memoryA_output;
    input wire [element_width * number_of_equations_per_cluster - 1 : 0]memoryP_output;
	input wire [element_width * no_of_units - 1 : 0] pKold_v2;
    input wire [element_width * no_of_units - 1 : 0]memoryR_output;
	input wire [memory_read_address_width-1:0]memoryR_read_address;
    input wire [number_of_equations_per_cluster * element_width - 1 : 0] memoryX_output;
	input wire memoryRprev_we;
	
	
	output wire [no_of_units * element_width - 1 : 0]memoryP_input;
    output wire [no_of_units * element_width - 1 : 0]memoryR_input;
    output wire [no_of_units * element_width - 1 : 0]memoryX_input;
    output wire finish;
	output wire mXv1_finish ;
      
	
	
	output wire result_mem_we_4 ; 
	output wire read_again;
	
	output wire result_mem_we_5;
	output wire[31:0] result_mem_counter_5;
	
	
	output wire result_mem_we_6;
	output wire read_again_2;
	
	output wire start;	 
	output wire vXv1_finish; 
	output wire finish_all;
	
	
	
    
    
    wire [element_width*no_of_units-1:0] rKold_prev;
   
    wire[31:0]rkold_read_address;
	
	
	
	//alu internal memories

    rKold_prev #(.no_of_units(no_of_units),.number_of_clusters(number_of_clusters),.number_of_equations_per_cluster(number_of_equations_per_cluster),.element_width (element_width ),.memories_address_width(memories_address_width)) 
	rk_prev(clk,memoryR_output,memoryR_read_address,rkold_read_address ,memoryRprev_we, rKold_prev);

        
	
	
	Alu #(.number_of_clusters(number_of_clusters),.number_of_equations_per_cluster(number_of_equations_per_cluster),.element_width (element_width ))
	alu(clk,reset,reset_vXv1,reset_mXv1,memoryA_output,memoryP_output,pKold_v2,memoryR_output,memoryX_output,rKold_prev,memoryP_input,memoryR_input,memoryX_input,finish,mXv1_finish,result_mem_we_4,rkold_read_address,result_mem_we_5,result_mem_counter_5,read_again,start,read_again_2,result_mem_we_6,vXv1_finish,finish_all);
	
endmodule