module memX(clk, input_data, write_enable, input_read_address, input_write_address, memory_output);
	
	parameter number_of_clusters = 20;
    parameter number_of_equations_per_cluster = 9;
    parameter element_width = 32;

	
	parameter no_of_units = 8;
	parameter memory_height = 1000;	
	parameter address_width= $clog2(memory_height)+1;
	
	input wire clk;
	input wire write_enable;
	input wire [no_of_units * element_width - 1 : 0] input_data;
	input wire [address_width - 1 : 0] input_write_address;
	input wire [address_width - 1 : 0] input_read_address;
	
	output wire [no_of_units * element_width - 1 : 0] memory_output;
	
	reg [no_of_units * element_width - 1 : 0] mem [0 : memory_height];
	// pragma attribute mem ram_block 1
	
	reg x=0;
	
	assign memory_output=mem[input_read_address];
	
	initial begin
		$readmemh("memx.txt", mem);
	end
	
	always @(posedge clk) 
		begin
			if( write_enable == 1'b1 ) 
				begin
					mem[input_write_address] <= input_data;
					x<=1;
				end
			
			 else if(x==1)
				 begin
					 $writememh("OUTPUTX.txt",mem);
					 x<=0;
				 end
			 end
			 
		 endmodule
