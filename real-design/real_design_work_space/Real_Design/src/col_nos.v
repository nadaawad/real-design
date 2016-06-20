module col_nos(clk,col_nos_read_address,col_nos_output);
	
	parameter no_of_elements_on_col_nos = 20 ;	 
	parameter no_of_row_by_vector_modules = 4; 
	
	input clk;
	input[31:0] col_nos_read_address;  // DONT MAKE THIS element_width
	output wire[no_of_row_by_vector_modules*no_of_elements_on_col_nos*32-1:0] col_nos_output; // DONT MAKE THIS element_width
	
	reg [no_of_row_by_vector_modules*no_of_elements_on_col_nos*32-1:0] mem [0 : 100000];	 // DONT MAKE THIS element_width
	
	
	assign col_nos_output = mem[col_nos_read_address];
		initial 
		begin
			$readmemh("col_nos.txt", mem);
		end		
		
		
	
	
endmodule	