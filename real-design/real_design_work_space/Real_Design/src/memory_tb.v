
`timescale 1 ns / 1 ps


module memory_tb (); 
	
	parameter no_of_elements_on_col_nos = 20 ;
	parameter no_of_elements_in_output = 8 ;
	parameter element_width = 32;   
	parameter no_of_units = 8;
	
	reg clk , read_preprocess,write_enable;	   
	reg [no_of_elements_on_col_nos*element_width-1:0] col_nos;	
	reg [31:0] no_of_multiples = 3;
	wire [no_of_elements_in_output*element_width-1:0] output_row; 
	
	

	P_Emap_8 P_E(clk,read_preprocess,write_enable,col_nos,output_row,no_of_multiples); 
	
	
	always begin
	#10
	clk <= ~ clk;
	end	
	
	initial
	begin
	  clk <=0; 
	  read_preprocess<=0;
	//  col_nos<= 256'h00000005000000530000005400000012000000110000001300000055FFFFFFFF;
	col_nos<=640'h0000601B00005FAE00005FAF000084F6000084F7000084740000609B0000601A0000609C000084F50000856300003B5D0000A9B30000AA350000AA340000A9470000A9480000A9B40000CE8FFFFFFFFF;
	  #50
	  read_preprocess<=1;
	  #20
	  read_preprocess<=0;
	  
	  #50
	  read_preprocess<=0;
	  
	end

endmodule
