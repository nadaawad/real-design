module multiples_fifo (clk,fifo_write_enable,fifo_read_address,fifo_write_address,fifo_input_data,fifo_output_data);

input clk,fifo_write_enable;

input [12:0]fifo_read_address;
input [12:0]fifo_write_address;
input[31:0] fifo_input_data;
output [31:0] fifo_output_data;

wire[12:0] real_write_address ;
assign real_write_address = (fifo_write_address !=0) ?(fifo_write_address-1) :13'h9 ;


reg [31:0] mem [0:100];


assign fifo_output_data = mem[fifo_read_address];

always @(posedge clk)
begin

	if(fifo_write_enable)
		begin
			mem[(fifo_write_address-1+10)%10] <= fifo_input_data;
		end

end





endmodule
