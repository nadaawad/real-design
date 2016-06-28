module TwoxOne_mux (mux_zero,mux_one,mux_output,mux_sel);

input wire [31:0] mux_zero;
input wire [31:0] mux_one;
input wire mux_sel;
output wire [31:0] mux_output;

assign mux_output = (mux_sel)?mux_one:mux_zero;


endmodule
