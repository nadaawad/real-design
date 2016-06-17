
module Eight_Organizer(clk,adder_row_input,start,adder_output);
parameter NI = 8;
input clk ;


input wire [NI*32-1:0] adder_row_input;
input wire start;

reg mux_select;
reg [31:0] mux_one = 32'b0;
wire [31:0] mux_output;
wire[31:0] pre_last_output;
output wire [31:0] adder_output ;

integer other_counter ;

TwoxOne_mux m1 (adder_output,mux_one,mux_output,mux_select);
EightxEight_Adder A1 (clk,adder_row_input,pre_last_output);
adder_subtractor final_adder (mux_output, pre_last_output, adder_output, 0,clk,1);


always @(posedge clk)
begin
	if(!start)
		begin
			other_counter<=0;	
			mux_select <= 1;
		end
else 
	begin
		other_counter <= other_counter +1 ;
		if(other_counter <9)
			begin
			mux_select <= 1;
			end
		else
			begin
				mux_select <= 0;
			end
	end
	
end






endmodule