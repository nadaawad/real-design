//-----------------------------------------------------------------------------
//
// Title       : row_by_vector
// Design      : improved cg
// Author      : Windows User
// Company     : nada
//
//-----------------------------------------------------------------------------
//
// File        : row by vector module.v
// Generated   : Wed Sep 16 23:46:23 2015
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
//{module {row_by_vector}}
module row_by_vector (clk,a,p,result);
	
//}} End of automatically maintained section

// -- Enter your statements here -- //

input clk;
wire clk;
input a,p;
wire[95:0] a,p;	
output result;
wire [31:0] result;	  



wire[31:0] m1_result,m2_result,m3_result,adder_1st_result,adder_2nd_result;



	  

multiply m1(a[95:64], p[95:64], clk, 1'b1,m1_result);

multiply m2 (a[63:32], p[63:32], clk, 1'b1, m2_result);

multiply m3 (a[31:0], p[31:0], clk, 1'b1, m3_result); 
adder_subtractor adder1 (m1_result,m2_result,adder_1st_result,1'b0,clk,1'b1); 
adder_subtractor adder2 (m3_result,32'b0,adder_2nd_result,1'b0,clk,1'b1);
adder_subtractor adder3 (adder_1st_result,adder_2nd_result,result,1'b0,clk,1'b1);




endmodule
