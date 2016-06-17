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
module adder_subtractor_with_control (A, B,result,op,clk,ce,outsider15,controlled_adder_output);
//}} End of automatically maintained section

// -- Enter your statements here -- //
input A,B;

wire[31:0] A,B;	  
wire [31:0]A_no_neg_zero,B_no_neg_zero;
assign A_no_neg_zero =(A[30:0]==31'b0)?32'b0:A;
assign B_no_neg_zero =(B[30:0]==31'b0)?32'b0:B;

input clk ,op,ce;
output result;	 

input wire outsider15;
output wire [31:0] controlled_adder_output;	
reg [31:0] previous_value;

wire [31:0] result;
wire[33:0] Radd,Rsub;
wire[33:0] modified_A , modified_B,modified_result ;

InputIEEE_8_23_to_8_23 A1 (.clk(clk),.rst(1'b0),.X(A_no_neg_zero), .R( modified_A));
InputIEEE_8_23_to_8_23 A2 (.clk(clk),.rst(1'b0),.X(B_no_neg_zero), .R( modified_B));

FPAddSub_8_23_uid2 adder(clk,~ce,modified_A ,modified_B, Radd, Rsub); 

assign modified_result= op? Rsub:Radd;

InputIEEE_8_23_to_8_23 A3 (.clk(clk),.rst(1'b0),.X(modified_result), .R( result));


assign controlled_adder_output = outsider15?result:previous_value;
always@(posedge clk)
	begin 
		if(outsider15==1)
			previous_value <= result;
		else if (!start )
			previous_value <= 32'b0;		
	end	

endmodule
