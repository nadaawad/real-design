module memR(clk, input_data, write_enable, input_read_address, input_write_address, memory_output,finish);
	

	parameter element_width = 32;
	parameter no_of_units = 8; 
	parameter memory_height = 1000;	
	parameter address_width= $clog2(memory_height)+1;
	
	input wire clk;
	input wire write_enable;
	input wire [no_of_units * element_width - 1 : 0] input_data;
	input wire [address_width - 1 : 0] input_write_address;
	input wire [address_width - 1 : 0] input_read_address; 
	
	
	output reg finish;
	output wire [no_of_units * element_width - 1 : 0] memory_output;
	
	reg x=0;
	reg [no_of_units * element_width - 1 : 0] mem [0 : memory_height];
	// pragma attribute mem ram_block 1
	
	assign memory_output=mem[input_read_address];
	
	initial 
		begin
			$readmemh("b.txt", mem);
			finish<=0;
		end
	
	always @(posedge clk) 
		begin
			if( write_enable == 1'b1 ) 
				begin
					mem[input_write_address] <= input_data;
				end
			end
			
	always@(posedge clk)
		 begin
			 if(x==1)
				 begin
					 x<=0;
				 end
			 end



endmodule
