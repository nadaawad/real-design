//-----------------------------------------------------------------------------
//
// Title       : vectorXvector_with_control
// Design      : real_design
// Author      : Windows User
// Company     : nada
//
//-----------------------------------------------------------------------------
//
// File        : vectorXvector_with_control.v
// Generated   : Fri Jun 17 14:19:32 2016
// From        : interface description file
// By          : Itf2Vhdl ver. 1.22
//
//-----------------------------------------------------------------------------
//
// Description : 
//
//-----------------------------------------------------------------------------
`timescale 1 ns / 1 ps


module vectorXvector_with_control (total,clk,reset,first_row_plus_additional,vector2,dot_product_output,finish,outsider_read_now); 

   
	parameter element_width=32;
	parameter no_of_units=8;
	

 
	integer counter2=0;

	
	
	input wire outsider_read_now;
	input wire clk;
	input wire reset;		
	input wire [element_width*no_of_units-1:0]vector2;
	input wire [element_width*no_of_units-1:0] first_row_plus_additional;
	input wire [31:0]total;
	
	output wire finish;		 
	
	output wire [element_width-1:0] dot_product_output ;
	 
	
	
	reg [no_of_units*element_width-1:0] first_row_input;
	reg [no_of_units*element_width-1:0] second_row_input;
	
	
	
	

	

	eight_Dot_Product_Multiply_with_control #(.no_of_units(no_of_units),.element_width(element_width))
	vXv(clk,reset ,first_row_input,second_row_input, dot_product_output,finish,outsider_read_now,total,I_am_ready);
		
				always @ (posedge clk)
					begin 
						if(reset)
							begin
							counter2<=0;
							
							end
						else if(!reset)
							begin  
								if(counter2 <total/no_of_units+2 && outsider_read_now)
									begin 
										
										first_row_input <= first_row_plus_additional;
										second_row_input <= vector2; 
										
										
										@(posedge clk);
										
										counter2 <=counter2+1;
									end	 
								
								end
							end


endmodule
