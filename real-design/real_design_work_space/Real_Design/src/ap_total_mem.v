module AP_total (clk, input_data,address,read_address,write_enable,memory_output);
	
	
	parameter element_width = 32;
	parameter no_of_units=8;
	
	
	input wire clk;
	input wire write_enable;
	input wire [element_width*no_of_units-1:0] input_data;
	input wire[31:0] address;
	input wire [31:0] read_address;
	

	output wire [element_width*no_of_units-1:0] memory_output;
	
	reg x=0;
	reg [element_width*no_of_units-1:0] mem [0 : 2000];
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