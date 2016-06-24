
module eight_Dot_Product_Multiply_with_control(total,clk,reset ,first_row_input,second_row_input, dot_product_output,finish,outsider_read_now );


parameter NI = 8;
integer repetition_times = (NI==8)?6:(NI==16)?6:0;

integer counter = 0;
integer ii=0; 
integer iii=0;
output reg finish ;

input wire [31:0]total;
input wire outsider_read_now;
input wire reset ;
input wire[32*NI-1:0] first_row_input;
input wire[32*NI-1:0] second_row_input;
reg save = 0;
reg adder_tree_start;
input clk ;
reg[32*NI-1:0] package_by_package;
wire [32*NI-1:0] multipliers_output_vector;
wire [31:0] adder_output;
output reg [31:0] dot_product_output;



wire[32*(NI/2)-1:0] demux_four_inputs;
reg demux_select;
reg flip; 
reg flip2;

//four_to_eight_demux demux_1(demux_four_inputs,demux_select,multipliers_output_vector);
N_to_2N_demux  #(.NI(NI))demux(demux_four_inputs,demux_select,multipliers_output_vector);
reg [32*(NI/2)-1:0] first_row_four_elements_subset;
reg [32*(NI/2)-1:0] second_row_four_elements_subset;	  

reg outsider1=0;
reg outsider2=0;
reg outsider3=0;
reg outsider4=0;  
reg outsider5=0;

wire outsider15;

genvar j ;
generate
for(j=0;j<NI/2;j=j+1) begin : instantiate_Multiplier

multiply m (first_row_four_elements_subset[32*(NI/2-j)-1-:32], second_row_four_elements_subset[32*(NI/2-j)-1-:32], clk, 1, demux_four_inputs[32*(NI/2-j)-1-:32]);
end
endgenerate


Eight_Organizer_with_control #(.NI(NI)) E_O (clk,package_by_package,adder_tree_start , adder_output,outsider4,outsider15);



always @ (posedge clk)
begin
	if(reset)
		begin
			
			finish <= 0;
			ii <=0;	   
			iii<=0;
		end
	else if(!reset) 
		begin
			if(ii < total/NI 	&& outsider5)
				begin
					package_by_package <= multipliers_output_vector;
					//@(posedge clk);
					ii <=ii+1;
				end
			else if(ii == total/NI)
				begin
					package_by_package <= 0; 
				end
			
		end
end							 

always @(posedge clk)
	begin

				if(iii < total/NI -1)
					begin 
						if(outsider15) 
							begin
								iii <= iii+1;  
							end	 
							
					end
				else if(iii == total/NI -1)
					begin 
						if(outsider15)
							begin	
								dot_product_output <= adder_output;
								finish<=1; 	  
							end
					end

	end




always @(posedge clk)
begin
if(reset) begin counter <= 0;
adder_tree_start <= 0;
end
else if(!reset)
begin
if(outsider_read_now) begin adder_tree_start <=1; end
counter <= counter+1;
end
end


always @(posedge clk)
	begin  
	outsider1 <= outsider_read_now;	
	outsider2 <= outsider1;
	outsider3 <= outsider2;
	outsider4 <= outsider3;
	outsider5<=  outsider4;
	 
	
	end	


always @(posedge clk)
	begin
		if(reset)
			begin
				//demux_select <= 0;
				flip <= 1;
			end
		else if(!reset && (outsider1 || ~flip))
			begin
				if(flip)
					begin
					//	demux_select <= 1;
						first_row_four_elements_subset <= first_row_input[NI*32-1-:(NI/2)*32];
						second_row_four_elements_subset <= second_row_input[NI*32-1-:(NI/2)*32];
						flip <= 0;
					end
				else if(!flip)
					begin
					//	demux_select <=0;
						first_row_four_elements_subset <= first_row_input[(NI/2)*32-1-:(NI/2)*32];
						second_row_four_elements_subset <= second_row_input[(NI/2)*32-1-:(NI/2)*32];
						flip	<= 1;
					end
			end
	end	
	
	always @(posedge clk)
	begin
		if(reset)
			begin
				demux_select <= 0;
				flip2 <= 1;
			end
		else if(!reset && (outsider3 || ~flip2))
			begin
				if(flip2)
					begin
						demux_select <= 1;
						flip2 <= 0;
					end
				else if(!flip2)
					begin
						demux_select <=0;
						flip2	<= 1;
					end
			end
	end


endmodule