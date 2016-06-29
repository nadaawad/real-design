module multiples_memory(clk,multiples_read_address,multiples_output); 
	parameter no_of_row_by_vector_modules = 4;
	parameter memory_A_height = 2000; 
	parameter address_width =$clog2(memory_A_height)+1;
	parameter multiples_memory_value_width=3;
	
	
	input wire clk;
	input wire[address_width-1:0] multiples_read_address;
	output wire [multiples_memory_value_width*no_of_row_by_vector_modules-1:0] multiples_output ; //dont make this element_width
	
	reg [multiples_memory_value_width*no_of_row_by_vector_modules-1 : 0] mem [0 : memory_A_height];	 // DONT MAKE THIS element_width
	
	
	assign multiples_output = mem[multiples_read_address];
		initial 
		begin
			$readmemh("multiples_matrix.txt", mem);
		end		
		
	

	
endmodule	