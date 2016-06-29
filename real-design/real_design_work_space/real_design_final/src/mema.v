module memA(clk,memA_read_address,memA_output,read_preprocess,no_of_multiples,I_am_ready);
	
	parameter no_of_elements_on_col_nos = 20 ;	 
	parameter no_of_row_by_vector_modules = 4; 	 
	parameter element_width =32;	 
	
    parameter no_of_elements_in_output = 8 ;	
	parameter overflow_allowed = no_of_elements_on_col_nos % no_of_elements_in_output;
	parameter additional =  (no_of_elements_in_output -  overflow_allowed)*element_width; 
	
	parameter memory_A_height = 2000;
	parameter address_width =$clog2(memory_A_height)+1;	 
	
	parameter multiples_memory_value_width =32;
	
	
	
	input wire read_preprocess;
	input wire[multiples_memory_value_width*no_of_row_by_vector_modules-1:0]no_of_multiples;
	input wire[no_of_row_by_vector_modules-1:0] I_am_ready ;
	
	reg[no_of_row_by_vector_modules-1:0] first_time_flag = -1;
	  
	

	output wire[no_of_row_by_vector_modules*element_width*no_of_elements_in_output-1:0] memA_output;
	
	input clk; 
	
	integer display_counter = 1;
	input[address_width-1:0] memA_read_address;  // DONT MAKE THIS element_width	
	
	wire[no_of_row_by_vector_modules*no_of_elements_on_col_nos*element_width-1:0] memA_extraction; // THIS MUST BE element_width
	
	reg[no_of_row_by_vector_modules*multiples_memory_value_width-1:0] fantastic_counters ;
	reg[no_of_row_by_vector_modules*multiples_memory_value_width-1:0] fantastic_counters_dash ;
	reg[no_of_row_by_vector_modules*multiples_memory_value_width-1:0] fantastic_counters_dash_dash ;
	

	wire [no_of_row_by_vector_modules*element_width*no_of_elements_in_output-1:0] value ;


	// THIS SHOULD BE AUTOMATED 
	
	assign value[4*element_width*no_of_elements_in_output-1-:element_width*no_of_elements_in_output]  = (no_of_elements_on_col_nos < fantastic_counters_dash_dash[4*multiples_memory_value_width-1-:multiples_memory_value_width]*no_of_elements_in_output)?
	({memA_extraction[(4*(no_of_elements_on_col_nos)-(fantastic_counters_dash_dash[4*multiples_memory_value_width-1-:multiples_memory_value_width]-1)*no_of_elements_in_output)*element_width-1-:overflow_allowed*element_width],{additional{1'b0}}})
	:memA_extraction[(4*(no_of_elements_on_col_nos)-(fantastic_counters_dash_dash[4*multiples_memory_value_width-1-:multiples_memory_value_width]-1)*no_of_elements_in_output)*element_width-1-:no_of_elements_in_output*element_width] 
	; 
	assign value[3*element_width*no_of_elements_in_output-1-:element_width*no_of_elements_in_output]  = (no_of_elements_on_col_nos < fantastic_counters_dash_dash[3*multiples_memory_value_width-1-:multiples_memory_value_width]*no_of_elements_in_output)?
	({memA_extraction[(3*(no_of_elements_on_col_nos)-(fantastic_counters_dash_dash[3*multiples_memory_value_width-1-:multiples_memory_value_width]-1)*no_of_elements_in_output)*element_width-1-:overflow_allowed*element_width],{additional{1'b0}}})
	:memA_extraction[(3*(no_of_elements_on_col_nos)-(fantastic_counters_dash_dash[3*multiples_memory_value_width-1-:multiples_memory_value_width]-1)*no_of_elements_in_output)*element_width-1-:no_of_elements_in_output*element_width] 
	;
	assign value[2*element_width*no_of_elements_in_output-1-:element_width*no_of_elements_in_output]  = (no_of_elements_on_col_nos < fantastic_counters_dash_dash[2*multiples_memory_value_width-1-:multiples_memory_value_width]*no_of_elements_in_output)?
	({memA_extraction[(2*(no_of_elements_on_col_nos)-(fantastic_counters_dash_dash[2*multiples_memory_value_width-1-:multiples_memory_value_width]-1)*no_of_elements_in_output)*element_width-1-:overflow_allowed*element_width],{additional{1'b0}}})
	:memA_extraction[(2*(no_of_elements_on_col_nos)-(fantastic_counters_dash_dash[2*multiples_memory_value_width-1-:multiples_memory_value_width]-1)*no_of_elements_in_output)*element_width-1-:no_of_elements_in_output*element_width] 
	;
	assign value[1*element_width*no_of_elements_in_output-1-:element_width*no_of_elements_in_output]  = (no_of_elements_on_col_nos < fantastic_counters_dash_dash[1*multiples_memory_value_width-1-:multiples_memory_value_width]*no_of_elements_in_output)?
	({memA_extraction[(1*(no_of_elements_on_col_nos)-(fantastic_counters_dash_dash[1*multiples_memory_value_width-1-:multiples_memory_value_width]-1)*no_of_elements_in_output)*element_width-1-:overflow_allowed*element_width],{additional{1'b0}}})
	:memA_extraction[(1*(no_of_elements_on_col_nos)-(fantastic_counters_dash_dash[1*multiples_memory_value_width-1-:multiples_memory_value_width]-1)*no_of_elements_in_output)*element_width-1-:no_of_elements_in_output*element_width] 
	;

	assign memA_output = value;
	
	reg [no_of_row_by_vector_modules*no_of_elements_on_col_nos*element_width-1:0] mem [0 : memory_A_height];	 // THIS MUST BE element_width
	
	
	assign memA_extraction = mem[memA_read_address];
		initial 
		begin
			$readmemh("A.txt", mem);	
			fantastic_counters[4*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			fantastic_counters[3*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			fantastic_counters[2*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			fantastic_counters[1*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			
			fantastic_counters_dash[4*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			fantastic_counters_dash[3*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			fantastic_counters_dash[2*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			fantastic_counters_dash[1*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			
			fantastic_counters_dash_dash[4*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			fantastic_counters_dash_dash[3*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			fantastic_counters_dash_dash[2*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
			fantastic_counters_dash_dash[1*multiples_memory_value_width-1-:multiples_memory_value_width] = 4'b0001;
		end		
		
	

		
	


			always @(posedge clk)
				begin
					if(read_preprocess || ~first_time_flag[4-1])
						begin
							@(posedge clk);	
							if(fantastic_counters[4*multiples_memory_value_width-1-:multiples_memory_value_width]<no_of_multiples[4*multiples_memory_value_width-1-:multiples_memory_value_width])
								begin
								//	@(I_am_ready[3]);
									fantastic_counters[4*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters[4*multiples_memory_value_width-1-:multiples_memory_value_width]+1;
									first_time_flag[4-1]<=0;
								end
							else 
								begin
									fantastic_counters[4*multiples_memory_value_width-1-:multiples_memory_value_width]<=1;
									first_time_flag[4-1]<=1;
								end

						end
				end
			always @(posedge clk)
				begin
					if(read_preprocess || ~first_time_flag[3-1])
						begin
							@(posedge clk);	
							if(fantastic_counters[3*multiples_memory_value_width-1-:multiples_memory_value_width]<no_of_multiples[3*multiples_memory_value_width-1-:multiples_memory_value_width])
								begin
								//	@(I_am_ready[2]);
									fantastic_counters[3*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters[3*multiples_memory_value_width-1-:multiples_memory_value_width]+1;
									first_time_flag[3-1]<=0;
								end
							else 
								begin
									fantastic_counters[3*multiples_memory_value_width-1-:multiples_memory_value_width]<=1;
									first_time_flag[3-1]<=1;
								end

						end
				end
			always @(posedge clk)
				begin
					if(read_preprocess || ~first_time_flag[2-1])
						begin
							@(posedge clk);	
							if(fantastic_counters[2*multiples_memory_value_width-1-:multiples_memory_value_width]<no_of_multiples[2*multiples_memory_value_width-1-:multiples_memory_value_width])
								begin
								//	@(I_am_ready[1]);
									fantastic_counters[2*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters[2*multiples_memory_value_width-1-:multiples_memory_value_width]+1;
									first_time_flag[2-1]<=0;
								end
							else 
								begin
									fantastic_counters[2*multiples_memory_value_width-1-:multiples_memory_value_width]<=1;
									first_time_flag[2-1]<=1;
								end

						end
				end
			always @(posedge clk)
				begin
					if(read_preprocess || ~first_time_flag[1-1])
						begin
							@(posedge clk);	
							if(fantastic_counters[1*multiples_memory_value_width-1-:multiples_memory_value_width]<no_of_multiples[1*multiples_memory_value_width-1-:multiples_memory_value_width])
								begin
								//	@(I_am_ready[0]);
									fantastic_counters[1*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters[1*multiples_memory_value_width-1-:multiples_memory_value_width]+1;
									first_time_flag[1-1]<=0;
								end
							else 
								begin
									fantastic_counters[1*multiples_memory_value_width-1-:multiples_memory_value_width]<=1;
									first_time_flag[1-1]<=1;
								end

						end

				end  
			
			always @(posedge clk)
				begin 
					fantastic_counters_dash[4*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters[4*multiples_memory_value_width-1-:multiples_memory_value_width];
					fantastic_counters_dash[3*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters[3*multiples_memory_value_width-1-:multiples_memory_value_width];
					fantastic_counters_dash[2*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters[2*multiples_memory_value_width-1-:multiples_memory_value_width];
					fantastic_counters_dash[1*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters[1*multiples_memory_value_width-1-:multiples_memory_value_width];
				end	
			always @(posedge clk)
				begin 
					fantastic_counters_dash_dash[4*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters_dash[4*multiples_memory_value_width-1-:multiples_memory_value_width];
					fantastic_counters_dash_dash[3*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters_dash[3*multiples_memory_value_width-1-:multiples_memory_value_width];
					fantastic_counters_dash_dash[2*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters_dash[2*multiples_memory_value_width-1-:multiples_memory_value_width];
					fantastic_counters_dash_dash[1*multiples_memory_value_width-1-:multiples_memory_value_width]<=fantastic_counters_dash[1*multiples_memory_value_width-1-:multiples_memory_value_width];
				end	
				

				
			
	
endmodule	