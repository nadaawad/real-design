
`timescale 1 ns / 1 ps


module tb_row_by_vector ();

		parameter element_width = 32;
	parameter memories_address_width=32;	
	
	parameter no_of_elements_on_col_nos = 20 ;	 
	parameter no_of_row_by_vector_modules = 4;	
	parameter no_of_units = no_of_row_by_vector_modules*2;
	parameter no_of_elements_in_p_emap_output = 8;	  
	
			reg	clk=0;	
			reg	reset =1; 
			reg	reset_mXv1 =0;	
			wire memories_pre_preprocess; 
			reg memories_preprocess=0;
	

	
	reg [memories_address_width - 1 : 0] memoryA_read_address;	
	reg [memories_address_width - 1 : 0] col_nos_read_address;	
	reg [memories_address_width - 1 : 0] multiples_read_address;	
	
	



	
	wire[no_of_row_by_vector_modules*element_width*no_of_units-1:0] memA_output;
	wire[no_of_row_by_vector_modules*no_of_elements_on_col_nos*32-1:0] col_nos_output;
	wire [32*no_of_row_by_vector_modules-1:0] multiples_output ;	
	wire [no_of_row_by_vector_modules*no_of_elements_in_p_emap_output*element_width-1:0] Emap_mem_output_row ;
	wire[31:0]total_with_additional_A;

	

   wire P_Emap_write_enable ;  
   wire [no_of_row_by_vector_modules-1:0] you_can_read;	
   wire[no_of_units*element_width-1:0] mXv1_result;


	
	genvar j;
		generate
		for(j=0;j<no_of_row_by_vector_modules;j=j+1) begin:instantiate_P_Emap_8	
			P_Emap_8 #(.no_of_units(no_of_units),.element_width(element_width),.no_of_elements_on_col_nos(no_of_elements_on_col_nos),.no_of_elements_in_output(no_of_elements_in_p_emap_output)) 
			P_Emap_mem (clk,memories_preprocess,P_Emap_write_enable,
			col_nos_output[(no_of_row_by_vector_modules-j)*no_of_elements_on_col_nos*32-1-:no_of_elements_on_col_nos*32],
			Emap_mem_output_row[((no_of_row_by_vector_modules-j))*no_of_elements_in_p_emap_output*element_width-1-:(no_of_elements_in_p_emap_output*element_width)],
			multiples_output[(no_of_row_by_vector_modules-j)*32-1-:32],you_can_read[no_of_row_by_vector_modules-j-1]);
		end
	endgenerate


	memA #(.no_of_elements_on_col_nos(no_of_elements_on_col_nos),.no_of_row_by_vector_modules(no_of_row_by_vector_modules),.element_width (element_width ))
	matA(clk,memoryA_read_address,memA_output,memories_preprocess,multiples_output);	 
	
	col_nos #(.no_of_elements_on_col_nos(no_of_elements_on_col_nos),.no_of_row_by_vector_modules(no_of_row_by_vector_modules))
	col_nos_memory(clk,col_nos_read_address,col_nos_output);
	
	multiples_memory #(.no_of_row_by_vector_modules(no_of_row_by_vector_modules))
	multiples_mat (clk,multiples_read_address,multiples_output);
	
	parameters_mem parameters_memory(total_with_additional_A);	
	
	
	matrix_by_vector_v3_with_control #(.no_of_row_by_vector_modules(no_of_units/2),.NI(no_of_units),.element_width (element_width ))
	mXv1_dash(clk,reset,reset_mXv1,memA_output,Emap_mem_output_row,mXv1_result,mXv1_finish,outsider_read_now,multiples_output,total_with_additional_A,memories_pre_preprocess,you_can_read);
	
	
	initial
		begin 
			clk<=0;	
			reset <=1; 
			reset_mXv1 <=0;	 
			
				 memoryA_read_address <=-1;
				 col_nos_read_address <=-1;
				 multiples_read_address<=-1;
			
			#60
			
			reset <=0; 
			reset_mXv1 <=1;	 
			
			#15

			
			#400

				 
			#190 
 
			
			#20 

			
			
			#2000
			$finish();
			
		end 
		
	always 
		begin 
			#10 clk <= ~ clk;
		end	 
		
	always @(posedge clk)
		begin
			if(memories_pre_preprocess)
				begin
				memoryA_read_address <= memoryA_read_address +1 ;
				 col_nos_read_address <= col_nos_read_address+1;
				 multiples_read_address<= multiples_read_address +1;
				end 	

		end	 
		
	always @(posedge clk)
		begin
		 memories_preprocess <= memories_pre_preprocess	;	
		 if(outsider_read_now)
			 begin
				 $display("%h",mXv1_result); 
			 end

		end	


endmodule
