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


module adder_subtractor_with_start (start,A, B,result,op,clk,ce,finish);

input A,B,start;
output reg finish=0;

wire[31:0] A,B;	 

wire [31:0]A_no_neg_zero,B_no_neg_zero;
assign A_no_neg_zero =(A[30:0]==31'b0)?32'b0:A;
assign B_no_neg_zero =(B[30:0]==31'b0)?32'b0:B;





input clk ,op,ce;
output result;
reg pip1=0;
reg pip2=0;
wire [31:0] result;


fpadd Adder(A_no_neg_zero, B_no_neg_zero, result,op,clk,1'b1);

always @(posedge clk)
	begin
			finish<=start;	
	end

endmodule
