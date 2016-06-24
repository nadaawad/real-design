module EightxEight_Adder_with_start(ExE_start,clk,inputs,summation,ExE_finish,ExE_finish_dash);
input wire clk;
parameter NI = 8;
input wire [NI*32-1:0] inputs ;
output wire[31:0] summation;   

input wire ExE_start;
reg ExE_start_dash=0;
reg ExE_start_dash_dash=0;

			   

wire[31:0] zero_one_sum ;
wire [31:0] two_three_sum;
wire [31:0] four_five_sum;
wire[31:0] six_seven_sum;
wire[31:0] zero_to_three_sum;
wire [31:0] four_to_seven_sum; 


wire[1:0] finish_first_two;
reg [1:0] finish_first_two_dash=0;
wire[1:0] finish_second_two;  
reg [1:0] finish_second_two_dash=0; 
wire[1:0] finish_third_two;  
reg [1:0] finish_third_two_dash=0;
output wire ExE_finish;	
output reg ExE_finish_dash;




//stage
 
adder_subtractor_with_start zero_one_adder (ExE_start_dash_dash,inputs[31:0], inputs[63:32], zero_one_sum, 1'b0,clk,1'b1,finish_first_two[1]);
adder_subtractor_with_start two_three_adder (ExE_start_dash_dash,inputs[95:64], inputs[127:96], two_three_sum, 1'b0,clk,1'b1,finish_first_two[0]);   

adder_subtractor_with_start four_five_adder (ExE_start_dash_dash,inputs[159:128], inputs[191:160], four_five_sum, 1'b0,clk,1'b1,finish_second_two[1]);
adder_subtractor_with_start six_seven_adder (ExE_start_dash_dash,inputs[223:192], inputs[255:224], six_seven_sum, 1'b0,clk,1'b1,finish_second_two[0]);

// stage
adder_subtractor_with_start zero_to_three_adder (&finish_first_two_dash,zero_one_sum, two_three_sum, zero_to_three_sum, 1'b0,clk,1'b1,finish_third_two[1]);
adder_subtractor_with_start four_to_seven_adder (&finish_second_two_dash,four_five_sum,six_seven_sum , four_to_seven_sum, 1'b0,clk,1'b1,finish_third_two[0]);
//stage
adder_subtractor_with_start zero_to_six_adder (&finish_third_two_dash,zero_to_three_sum, four_to_seven_sum, summation, 1'b0,clk,1'b1,ExE_finish);


always@(posedge clk)
	begin
		if(finish_first_two!=0 && ! (&finish_first_two_dash))
			begin
				finish_first_two_dash <= finish_first_two_dash | finish_first_two;	
			end	
		else if(&finish_first_two_dash && !(&finish_first_two))	
			begin 
				finish_first_two_dash<=0;
			end	
	end
	
always@(posedge clk)
	begin
		if(finish_second_two!=0 && ! (&finish_second_two_dash))
			begin
				finish_second_two_dash <= finish_second_two_dash | finish_second_two;	
			end	
		else if(&finish_second_two_dash && !(&finish_second_two))	
			begin 
				finish_second_two_dash<=0;
			end	
	end	
	
always@(posedge clk)
	begin
		if(finish_third_two!=0 && ! (&finish_third_two_dash))
			begin
				finish_third_two_dash <= finish_third_two_dash | finish_third_two;	
			end	
		else if(&finish_third_two_dash && !(&finish_third_two))	
			begin 
				finish_third_two_dash<=0;
			end	
	end
	
always @(posedge clk)
	begin
		ExE_finish_dash<=ExE_finish;
		
		ExE_start_dash<=ExE_start;
		ExE_start_dash_dash <=ExE_start_dash;
	end	


endmodule
