module multiples_memory(clk,multiples_read_address,multiples_output); 
	parameter no_of_row_by_vector_modules = 4;
	
	input wire clk;
	input wire[31:0] multiples_read_address;
	output wire [32*no_of_row_by_vector_modules-1:0] multiples_output ; //dont make this element_width
	
	reg [32*no_of_row_by_vector_modules-1 : 0] mem [0 : 100000];	 // DONT MAKE THIS element_width
	
	
	assign multiples_output = mem[multiples_read_address];
		initial 
		begin
			$readmemh("multiples_matrix.txt", mem);
		end		
		
	

	
endmodule	