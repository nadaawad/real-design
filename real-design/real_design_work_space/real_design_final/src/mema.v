module memA(clk,memA_read_address,memA_output,read_preprocess,no_of_multiples,I_am_ready);
	
	parameter no_of_elements_on_col_nos = 20 ;	 
	parameter no_of_row_by_vector_modules = 4; 	 
	parameter element_width =32;	 
	parameter no_of_units = no_of_row_by_vector_modules * 2;  
	parameter overflow_allowed = no_of_elements_on_col_nos % no_of_units;
	parameter additional =  (no_of_units -  overflow_allowed)*element_width; 
	
	parameter memory_A_height = 2000;
	parameter address_width =$clog2(memory_A_height)+1;	 
	
	parameter multiples_memory_value_width =32;
	
	input wire read_preprocess;
	input wire[multiples_memory_value_width*no_of_row_by_vector_modules-1:0]no_of_multiples;
	input wire[no_of_row_by_vector_modules-1:0] I_am_ready ;
	
	  
	reg first_time_flag_4 =1;
	reg first_time_flag_3 =1;
	reg first_time_flag_2 =1;
	reg first_time_flag_1 =1;
		
	output wire[no_of_row_by_vector_modules*element_width*no_of_units-1:0] memA_output;
	
	input clk; 
	
	integer display_counter = 1;
	input[address_width-1:0] memA_read_address;  // DONT MAKE THIS element_width	
	
	wire[no_of_row_by_vector_modules*no_of_elements_on_col_nos*element_width-1:0] memA_extraction; // THIS MUST BE element_width
	
	integer i=1; 
	integer j=1;
	integer k=1;
	integer q=1; 
	
	integer ii;
	integer jj;
	integer kk;
	integer qq;	
	
	integer iii;
	integer jjj;
	integer kkk;
	integer qqq;  
	
	wire[element_width*no_of_units-1:0] value_4; 
	wire[element_width*no_of_units-1:0] value_3;
	wire[element_width*no_of_units-1:0] value_2;
	wire[element_width*no_of_units-1:0] value_1;

	// THIS SHOULD BE AUTOMATED 
	
	assign value_4  = (no_of_elements_on_col_nos < iii*no_of_units)?
	({memA_extraction[(4*(no_of_elements_on_col_nos)-(iii-1)*no_of_units)*element_width-1-:overflow_allowed*element_width],{additional{1'b0}}})
	:memA_extraction[(4*(no_of_elements_on_col_nos)-(iii-1)*no_of_units)*element_width-1-:no_of_units*element_width] 
	; 
	assign value_3  = (no_of_elements_on_col_nos < jjj*no_of_units)?
	({memA_extraction[(3*(no_of_elements_on_col_nos)-(jjj-1)*no_of_units)*element_width-1-:overflow_allowed*element_width],{additional{1'b0}}})
	:memA_extraction[(3*(no_of_elements_on_col_nos)-(jjj-1)*no_of_units)*element_width-1-:no_of_units*element_width] 
	;
	assign value_2  = (no_of_elements_on_col_nos < kkk*no_of_units)?
	({memA_extraction[(2*(no_of_elements_on_col_nos)-(kkk-1)*no_of_units)*element_width-1-:overflow_allowed*element_width],{additional{1'b0}}})
	:memA_extraction[(2*(no_of_elements_on_col_nos)-(kkk-1)*no_of_units)*element_width-1-:no_of_units*element_width] 
	;
	assign value_1  = (no_of_elements_on_col_nos < qqq*no_of_units)?
	({memA_extraction[(1*(no_of_elements_on_col_nos)-(qqq-1)*no_of_units)*element_width-1-:overflow_allowed*element_width],{additional{1'b0}}})
	:memA_extraction[(1*(no_of_elements_on_col_nos)-(qqq-1)*no_of_units)*element_width-1-:no_of_units*element_width] 
	;

	assign memA_output = {value_4,value_3,value_2,value_1};
	
	reg [no_of_row_by_vector_modules*no_of_elements_on_col_nos*element_width-1:0] mem [0 : memory_A_height];	 // THIS MUST BE element_width
	
	
	assign memA_extraction = mem[memA_read_address];
		initial 
		begin
			$readmemh("A.txt", mem);
		end		
		
	

		
	


			always @(posedge clk)
				begin
					if(read_preprocess || ~first_time_flag_4)
						begin
							@(posedge clk);	
							if(i<no_of_multiples[4*32-1-:32])
								begin
								//	@(I_am_ready[3]);
									i<=i+1;
									first_time_flag_4<=0;
								end
							else 
								begin
									i<=1;
									first_time_flag_4<=1;
								end

						end
				end
			always @(posedge clk)
				begin
					if(read_preprocess || ~first_time_flag_3)
						begin
							@(posedge clk);	
							if(j<no_of_multiples[3*32-1-:32])
								begin
								//	@(I_am_ready[2]);
									j<=j+1;
									first_time_flag_3<=0;
								end
							else 
								begin
									j<=1;
									first_time_flag_3<=1;
								end

						end
				end
			always @(posedge clk)
				begin
					if(read_preprocess || ~first_time_flag_2)
						begin
							@(posedge clk);	
							if(k<no_of_multiples[2*32-1-:32])
								begin
								//	@(I_am_ready[1]);
									k<=k+1;
									first_time_flag_2<=0;
								end
							else 
								begin
									k<=1;
									first_time_flag_2<=1;
								end

						end
				end
			always @(posedge clk)
				begin
					if(read_preprocess || ~first_time_flag_1)
						begin
							@(posedge clk);	
							if(q<no_of_multiples[1*32-1-:32])
								begin
								//	@(I_am_ready[0]);
									q<=q+1;
									first_time_flag_1<=0;
								end
							else 
								begin
									q<=1;
									first_time_flag_1<=1;
								end

						end

				end  
			
			always @(posedge clk)
				begin 
					ii<=i;
					jj<=j;
					kk<=k;
					qq<=q;
				end	
			always @(posedge clk)
				begin 
					iii<=ii;
					jjj<=jj;
					kkk<=kk;
					qqq<=qq;
				end	
				

				
		/*	always @(posedge clk)
				begin
				$display(" DISPLAY COUNTER %d",display_counter);
				$display(" Read_Preprocess %b",read_preprocess);
			

				$display("%h",value_4);
				$display("%h",value_3);
				$display("%h",value_2);
				$display("%h",value_1);
	
			   	display_counter =  display_counter+1;	  
			
				end		
		*/
			
	
endmodule	