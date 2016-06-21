
`timescale 1 ns / 1 ps


module row_by_vector_with_control (clk,a,p,result,give_me_only,no_of_multiples,start_row_by_vector,decoder_read_now,reset,you_can_read);



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
reg pipeline0=0,pipeline1=0,pipeline2=0,pipeline3=0,pipeline4=0,pipeline5=0,pipeline6=0;	  
input wire reset;
integer i = 0;	   

reg initialization_counter = 0;		
input wire you_can_read;





	  

//multiply m1(a[95:64], p[95:64], clk, 1'b1,m1_result);
//multiply m2 (a[63:32], p[63:32], clk, 1'b1, m2_result);
//multiply m3 (a[31:0], p[31:0], clk, 1'b1, m3_result); 
//adder_subtractor adder1 (m1_result,m2_result,adder_1st_result,1'b0,clk,1'b1); 
//adder_subtractor adder2 (m3_result,32'b0,adder_2nd_result,1'b0,clk,1'b1);
//adder_subtractor adder3 (adder_1st_result,adder_2nd_result,result,1'b0,clk,1'b1);	

eight_Dot_Product_Multiply_with_control_row  edomwcr(clk,start_row_by_vector ,a,p, result,dot_product_finish,you_can_read,no_of_multiples,prepare_my_new_input);

// note if you don't need to accelerate the special case where #of multiples =1 , you don't need the reset signal
	// nor the total and no_of_units parameters. 
	
// note , both the special and slower case seem to work well now , although I am more relieved with the non-special case
	// also for the non special case , decoder_read_now<=pipeline5 will work right , but won't work for the special case
// overall it seems logical , and i didn't do any fabraka , i guess .		

always @(posedge clk) 
	begin 
	if (reset) begin  i<=0 ;end
	else if(prepare_my_new_input)
		begin  
		   give_me_only<=1;	
		end	
	else give_me_only<=0;	
		
	end	  
	
	always @(posedge clk)
		begin
			if (prepare_my_new_input)	   
			begin 
				pipeline0 <=1 ;
			end 
			else 
				begin 
					pipeline0 <= 0 	 ;
				end	
			
			pipeline1 <=pipeline0;
			//pipeline2 <=pipeline1;
			//pipeline3 <= pipeline2;
			//pipeline4 <= pipeline3;
			//pipeline5 <=pipeline4;
			//decoder_read_now <= pipeline5;
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
