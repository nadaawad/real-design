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
module decoder (in ,sel,out);
	parameter no_of_units =4;
	parameter element_width =32;
	input in;
	wire[element_width*no_of_units-1:0] in;
	input sel;
	output out;
	reg[2*element_width*no_of_units-1:0] out;
	
	always @(in or sel)
		begin
			if(!sel)
				out[element_width*no_of_units-1:0]<=in;
			else
				out[2*element_width*no_of_units-1:element_width*no_of_units] <=in;
				
			end
			
	
	
endmodule
