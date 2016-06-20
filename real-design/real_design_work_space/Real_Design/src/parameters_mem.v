module parameters_mem(total_with_additional_A);
	
	
	output wire[31:0] total_with_additional_A ;
	reg [31:0]mem[0:10];
	assign total_with_additional_A= mem[0];
	
	initial
		begin
			$readmemh("Parameters.txt", mem);	
		end	
	
	
endmodule	