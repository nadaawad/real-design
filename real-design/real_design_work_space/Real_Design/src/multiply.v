//-----------------------------------------------------------------------------
//
// Title       : multipy
// Design      : cluster jacoubi
// Author      : Windows User
// Company     : nada
//
//-----------------------------------------------------------------------------
//
// File        : mutiply.v
// Generated   : Fri Sep 18 15:01:08 2015
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
//{module {multipy}}
module multiply (A, B, clk, ce, result);
//}} End of automatically maintained section

// -- Enter your statements here -- //
input A,B;
wire[31:0] A,B;
input clk,ce;
output result;
wire [31:0] result;
wire[33:0] modified_A , modified_B ,out;

InputIEEE_8_23_to_8_23 A1 (.clk(clk),.rst(1'b0),.X(A), .R( modified_A));
InputIEEE_8_23_to_8_23 A2 (.clk(clk),.rst(1'b0),.X(B), .R( modified_B));

FPMultiplier_8_23_8_23_8_23_uid2 m(.clk(clk),.rst(~ce),.X(modified_A),.Y(modified_B),.R(out)); 

OutputIEEE_8_23_to_8_23 Out1 (.clk(clk), .rst(1'b0),.X(out),.R(result));


endmodule
