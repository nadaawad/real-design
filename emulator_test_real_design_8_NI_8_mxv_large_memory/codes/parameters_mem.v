module parameters_mem(total_with_additional_A,total);
	
	
	output wire[31:0] total_with_additional_A ;
	output wire[31:0] total;
	reg [31:0]mem[0:1];
	assign total_with_additional_A= mem[0];
	assign total= mem[1];
	
	initial
		begin
			$readmemh("Parameters.txt", mem);	
		end	
	
	
endmodule	