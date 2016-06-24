module multiples_fifo (clk,fifo_write_enable,fifo_read_address,fifo_write_address,fifo_input_data,fifo_output_data);

input clk,fifo_write_enable;

input [3:0]fifo_read_address;
input [3:0]fifo_write_address;
input[31:0] fifo_input_data;
output [31:0] fifo_output_data;


reg [31:0] mem [1:10];


assign fifo_output_data = mem[fifo_read_address];

always @(posedge clk)
begin

	if(fifo_write_enable)
		begin
			mem[fifo_write_address] <= fifo_input_data;
		end

end





endmodule
