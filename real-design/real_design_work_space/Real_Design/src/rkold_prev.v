module rKold_prev(clk, input_data, address,read_address ,write_enable,memory_output);
	
	parameter number_of_clusters = 1;
	parameter number_of_equations_per_cluster = 9;
	parameter element_width = 32;
	parameter address_width = 20;
	parameter memories_address_width=20;
	parameter no_of_units = 8; 
	
	input wire clk;
	input wire write_enable;
	input wire [element_width*no_of_units-1:0] input_data;
	input wire [address_width-1:0] address; 
	input wire [address_width-1:0] read_address;
	

	output wire [element_width*no_of_units-1:0] memory_output;
	
	reg x=0;
	reg [element_width*no_of_units-1:0] mem [0 : 1000];
	// pragma attribute mem ram_block 1
	
	assign memory_output=mem[read_address];
	
	
	
	always @(posedge clk) 
		begin
			if( write_enable == 1'b1 ) 
				begin
					mem[address] <= input_data; 
				end

			end
			
	

endmodule