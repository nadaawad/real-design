//-----------------------------------------------------------------------------
//
// Title       : adder_subtractor
// Design      : cluster jacoubi
// Author      : Windows User
// Company     : nada
//
//-----------------------------------------------------------------------------
//
// File        : adder subtractor.v
// Generated   : Fri Sep 18 15:00:44 2015
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
//{module {adder_subtractor}}
module adder_subtractor (A, B,result,op,clk,ce);
//}} End of automatically maintained section

// -- Enter your statements here -- //
input A,B;

wire[31:0] A,B;	 

wire [31:0]A_no_neg_zero,B_no_neg_zero;
assign A_no_neg_zero =(A[30:0]==31'b0)?32'b0:A;
assign B_no_neg_zero =(B[30:0]==31'b0)?32'b0:B;

input clk ,op,ce;
output result;

wire [31:0] result;

fpadd Adder(A_no_neg_zero, B_no_neg_zero, result,op,clk,1'b1);

endmodule
