module EightxEight_Adder(clk,inputs,summation);
input wire clk;
parameter NI = 8;
input wire [NI*32-1:0] inputs ;
output wire[31:0] summation;


wire[31:0] zero_one_sum ;
wire [31:0] two_three_sum;
wire [31:0] four_five_sum;
wire[31:0] six_seven_sum;
wire[31:0] zero_to_three_sum;
wire [31:0] four_to_seven_sum;

//stage
 
adder_subtractor zero_one_adder (inputs[31:0], inputs[63:32], zero_one_sum, 1'b0,clk,1'b1);
adder_subtractor two_three_adder (inputs[95:64], inputs[127:96], two_three_sum, 1'b0,clk,1'b1);
adder_subtractor four_five_adder (inputs[159:128], inputs[191:160], four_five_sum, 1'b0,clk,1'b1);
adder_subtractor six_seven_adder (inputs[223:192], inputs[255:224], six_seven_sum, 1'b0,clk,1'b1);

// stage
adder_subtractor zero_to_three_adder (zero_one_sum, two_three_sum, zero_to_three_sum, 1'b0,clk,1'b1);
adder_subtractor four_to_seven_adder (four_five_sum,six_seven_sum , four_to_seven_sum, 1'b0,clk,1'b1);
//stage
adder_subtractor zero_to_six_adder (zero_to_three_sum, four_to_seven_sum, summation, 1'b0,clk,1'b1);


endmodule
