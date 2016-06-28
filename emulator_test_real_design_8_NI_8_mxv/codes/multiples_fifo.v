module multiples_fifo (clk,fifo_write_enable,fifo_read_address,fifo_write_address,fifo_input_data,fifo_output_data);

parameter multiples_memory_value_width = 3;
input clk,fifo_write_enable;

input [3:0]fifo_read_address;
input [3:0]fifo_write_address;
input[multiples_memory_value_width-1:0] fifo_input_data;
output [multiples_memory_value_width-1:0] fifo_output_data;

wire[3:0] real_write_address ;
assign real_write_address = (fifo_write_address !=0) ?(fifo_write_address-1) :4'h9 ;


reg [multiples_memory_value_width-1:0] mem [0:10];


assign fifo_output_data = mem[fifo_read_address];

always @(posedge clk)
begin

	if(fifo_write_enable)
		begin
			mem[(fifo_write_address-1+10)%10] <= fifo_input_data;
		end

end





endmodule
