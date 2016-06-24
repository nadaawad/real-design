module N_to_2N_demux (in ,sel,out);
	parameter no_of_units=8;
	parameter element_width =32;
	input in;
	wire[element_width*(no_of_units/2)-1:0] in;
	input sel;
	output out;
	reg[2*element_width*(no_of_units/2)-1:0] out;
	
	always @(in or sel)
		begin
			if(!sel)
				out[element_width*(no_of_units/2)-1:0]<=in;
			else
				out[2*element_width*(no_of_units/2)-1:element_width*(no_of_units/2)] <=in;
				
			end
	
endmodule
