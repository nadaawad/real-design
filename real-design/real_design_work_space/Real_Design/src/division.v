//-----------------------------------------------------------------------------
//
// Title       : division
// Design      : basic conjugate gradient
// Author      : Windows User
// Company     : nada
//
//-----------------------------------------------------------------------------
//
// File        : division.v
// Generated   : Thu Aug 27 15:02:59 2015
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
//{module {division}}
module division ( clk ,start ,A ,B ,result ,finish_flag );

output [31:0] result ;
wire [31:0] result ;
output finish_flag ;
reg  finish_flag ;

input clk ;
wire clk ;
input start ;
wire start ;
input [31:0] A ;
wire [31:0] A ;
input [31:0] B ;
wire [31:0] B ;
reg count;




FP_div D(A, B, clk, start, result);


initial
	begin
		count=0;
		finish_flag=0;
		
	end	
	
	always@(posedge clk)
		begin 
			if(start==0)
				begin
					count<=0;
	                finish_flag<=0;
		         end
		
				
			else if(start&&count<1'b1) 
				begin
					count<= count+1;
					
				end
				
			else if(start&&count==1'b1)
				begin
					finish_flag <= 1; 
				
				end
				
				
			end	
//}} End of automatically maintained section

// -- Enter your statements here -- //

endmodule
