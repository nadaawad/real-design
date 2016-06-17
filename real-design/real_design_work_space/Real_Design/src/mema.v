module memA( input_address,clk, memory_output);

parameter number_of_clusters = 1;
parameter number_of_equations_per_cluster = 9;
parameter element_width = 32;
parameter address_width = 20;
parameter memories_address_width=20;

input wire clk;
input wire [address_width - 1 : 0] input_address;

output wire [element_width*(3* number_of_equations_per_cluster-2*2+2)-1 : 0] memory_output;

reg [element_width*(3* number_of_equations_per_cluster-2*2+2)-1 : 0] mem [0 : number_of_clusters - 1];
// pragma attribute mem ram_block 1

assign memory_output=mem[input_address];

initial 
	begin
		$readmemh("A.txt", mem);
	end
  
endmodule
