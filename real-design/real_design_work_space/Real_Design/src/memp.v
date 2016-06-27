module memP(clk, input_data, write_enable, input_read_address, input_write_address, memory_output,finish);
	
	parameter number_of_clusters = 1;
	parameter number_of_equations_per_cluster = 9;
	parameter element_width = 32;
	parameter address_width = 20;
	parameter memories_address_width=20;
	
	input wire clk;
	input wire write_enable;
	input wire [number_of_equations_per_cluster * element_width - 1 : 0] input_data;
	input wire [address_width - 1 : 0] input_write_address;
	input wire [address_width - 1 : 0] input_read_address;
	
	output reg finish;
	output wire [number_of_equations_per_cluster * element_width - 1 : 0] memory_output;
	
	reg x=0;
	reg [number_of_equations_per_cluster * element_width - 1 : 0] mem [0 : number_of_clusters - 1];
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
					 //$writememh("OUTPUTX.txt",mem);
					 x<=0;
				 end
			 end

endmodule
