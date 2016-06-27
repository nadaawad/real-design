
`timescale 1 ns / 1 ps


module row_by_vector_with_control (clk,a,p,result,give_me_only,no_of_multiples,start_row_by_vector,decoder_read_now,reset,you_can_read,I_am_ready);


parameter NI=8;
parameter element_width = 32;


input clk;
wire clk;
input a,p;
wire[NI*element_width-1:0] a,p;	
output result;
wire [element_width-1:0] result;	  


input wire [31:0] no_of_multiples;
integer no_of_multiples_counter = 0;
reg counter_for_deasserting=0;

//wire[element_width-1:0] m1_result,m2_result,m3_result,adder_1st_result,adder_2nd_result;	 

	
output reg give_me_only =0;  // this rises to one when ALL multiples have been calculated
// for example if we have 8 multipliers and the row has 24 entries , this rises after all the three multiples are done
output reg decoder_read_now;		
input wire start_row_by_vector;
reg pipeline0=0,pipeline1=0;	  
input wire reset;
integer i = 0;	   

reg initialization_counter = 0;		
input wire you_can_read;
wire fake_prepare0;

output wire I_am_ready;




eight_Dot_Product_Multiply_with_control_row  edomwcr(clk,reset,start_row_by_vector ,a,p, result,dot_product_finish,you_can_read,no_of_multiples,prepare_my_new_input,fake_prepare0,I_am_ready);

	

always @(posedge clk) 
	begin 
	if(!reset && prepare_my_new_input)
		begin  
		   give_me_only<=1;	
		end	
	else give_me_only<=0;	
		
	end	  
	
	always @(posedge clk)
		begin
			if (fake_prepare0)	   
			begin 
				pipeline0 <=1 ;
			end 
			else 
				begin 
					pipeline0 <= 0 	 ;
				end	
			
			pipeline1 <=pipeline0;
			decoder_read_now <= pipeline1;
			
		end	 
		
		always @(posedge clk)
			begin 
				if(reset)
					initialization_counter<=0;
				else if(start_row_by_vector)
					initialization_counter<=1;
				
			end 	
	

		

endmodule
