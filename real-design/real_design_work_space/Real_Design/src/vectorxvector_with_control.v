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


module vectorXvector_with_control (clk,reset,first_row_plus_additional,vector2,result,finish,outsider_read_now); 

    parameter number_of_equations_per_cluster=16;
	parameter element_width_modified=34;
	parameter element_width=32;
	parameter no_of_units=8;
	parameter count=number_of_equations_per_cluster/no_of_units;
	parameter additional = no_of_units-(number_of_equations_per_cluster%no_of_units); 
	parameter total = number_of_equations_per_cluster+additional ;
	parameter number_of_clusters=1;
	integer counter = 0;
	integer i=0; 
	integer counter2=0;
	integer counter3=0;
	
	
	input wire outsider_read_now;
	input wire clk;
	input wire reset;		
	input wire [element_width*no_of_units-1:0]vector2;
	input wire [element_width*no_of_units-1:0] first_row_plus_additional;
	//input wire [element_width*total-1:0] AP_total ;
	
	
	
	//output wire [element_width*number_of_equations_per_cluster-1:0] AP;
	output wire [element_width-1:0]result;
	output wire finish;		 
	//output reg vector1_mem_we;
	 
	
	
	reg [no_of_units*element_width-1:0] first_row_input;
	reg [no_of_units*element_width-1:0] second_row_input;
	
	
	
	
	//assign AP=AP_total[element_width*total-1:element_width*total-element_width*number_of_equations_per_cluster];
	

	eight_Dot_Product_Multiply_with_control #(.NOE(number_of_equations_per_cluster))
	vXv(clk,reset,first_row_input,second_row_input, result,finish,outsider_read_now);
		
//fiveonetwo_Dot_Product_Multiply#(.NOE(number_of_equations_per_cluster))
//vXv(clk,reset,first_row_input,second_row_input, result,finish );

//Sixteen_Dot_Product_Multiply #(.NOE(number_of_equations_per_cluster))
//vXv(clk,reset,first_row_input,second_row_input, result,finish );

//onezerotwofour_Dot_Product_Multiply #(.NOE(number_of_equations_per_cluster))
	//vXv(clk,reset,first_row_input,second_row_input, result,finish );
	
	initial
		begin
		
		  counter2<=0;
		end
		



	always @(posedge clk)
		begin
			if(reset)
				begin
				counter<=0;
				
				end
			else if(!reset)
				begin
					if (counter==0)
						begin 
							
							
						end
						counter <= counter+1;
					end
				end
				
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
