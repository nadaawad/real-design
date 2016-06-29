//-----------------------------------------------------------------------------
//
// Title       : decoder
// Design      : cluster jacoubi
// Author      : Windows User
// Company     : nada
//
//-----------------------------------------------------------------------------
//
// File        : decoder.v
// Generated   : Fri Oct  2 15:17:38 2015
// From        : interface description file
// By          : Itf2Vhdl ver. 1.22
//
//-----------------------------------------------------------------------------
//
// Description : 
//
//-----------------------------------------------------------------------------
`timescale 1 ns / 1 ps

//{{ Section below this comment is automatically maintained
//   and may be overwritten
//{module {decoder}}
module decoder_with_control (clk,in ,out,decoder_read_now,outsider_read_now);
	parameter no_of_row_by_vector_modules =4;
	parameter element_width =32;
	input in; 
	wire[element_width*no_of_row_by_vector_modules-1:0] in;
	input wire clk ;
	output out;	
	input wire decoder_read_now;
	reg[2*element_width*no_of_row_by_vector_modules-1:0] out;  
	reg internal_select;  
	
	output reg outsider_read_now;
	reg delay_flag=0;							   
	initial 
		begin
		internal_select <=1;
		outsider_read_now<=0;
		end
	
	always @(posedge clk)
		begin 
			if(decoder_read_now)
				begin 
					if(!internal_select)
						begin
								out[element_width*no_of_row_by_vector_modules-1:0]<=in;	
								internal_select <=1;
							//	delay_flag<=1;	  
							outsider_read_now <=1;
								
						end
					else 
						begin 
							out[2*element_width*no_of_row_by_vector_modules-1:element_width*no_of_row_by_vector_modules] <=in;	
							internal_select <=0;
							//delay_flag <=0;  
							outsider_read_now <=0;
						end	  
				end	
			else begin 
					outsider_read_now<=0;
				 end
		end
			
	
	
endmodule
