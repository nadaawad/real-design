
module eight_Dot_Product_Multiply_with_control(clk,reset ,first_row_input,second_row_input, dot_product_output,finish,outsider_read_now,total,I_am_ready);


parameter element_width=32;
parameter no_of_units = 8;
 

reg initialization_counter = 1;
reg initialization_counter2 = 1;  


integer ii=0; 
integer iii=0; 

output reg finish ;	 



input wire [31:0] total; // THIS IS Equivalent to the number of multiples in the  "edpmwc_row" module . but this has to be devided by no_of_units


input wire outsider_read_now;
input wire reset ;
input wire[element_width*no_of_units-1:0] first_row_input;
input wire[element_width*no_of_units-1:0] second_row_input;

reg [element_width*no_of_units-1:0] anding_mask = {no_of_units{32'h80000000}}; // IF 64 CHANGE MANUALY  
output reg I_am_ready;


reg save = 0;
reg adder_tree_start=0;
input clk ;
reg[element_width*no_of_units-1:0] package_by_package;
wire [element_width*no_of_units-1:0] multipliers_output_vector;
wire [element_width-1:0] adder_output;
output reg [element_width-1:0] dot_product_output;


wire[element_width*(no_of_units/2)-1:0] demux_four_inputs;
reg demux_select;
reg flip=1; 
reg flip2=1;

//four_to_eight_demux demux_1(demux_four_inputs,demux_select,multipliers_output_vector);
N_to_2N_demux  #(.no_of_units(no_of_units))demux(demux_four_inputs,demux_select,multipliers_output_vector);
reg [element_width*(no_of_units/2)-1:0] first_row_four_elements_subset;
reg [element_width*(no_of_units/2)-1:0] second_row_four_elements_subset;	  

reg outsider1=0;
reg outsider2=0;
reg outsider3=0;
reg outsider4=0;  
reg outsider5=0;   




wire ExE_finish;
wire final_adder_finish_dash;

genvar j ;
generate
for(j=0;j<no_of_units/2;j=j+1) begin : instantiate_Multiplier

multiply m (first_row_four_elements_subset[element_width*(no_of_units/2-j)-1-:element_width], second_row_four_elements_subset[element_width*(no_of_units/2-j)-1-:element_width], clk, 1, demux_four_inputs[element_width*(no_of_units/2-j)-1-:element_width]);
end
endgenerate


Eight_Organizer_with_control_row #(.no_of_units(no_of_units)) E_O (clk,package_by_package,adder_tree_start , adder_output,outsider4,final_adder_finish_dash,ExE_finish);


always @(negedge clk)
	begin
		if(outsider1)
			begin
							
				if  (
				(((first_row_input & anding_mask)^(second_row_input & anding_mask))== anding_mask) 
									||
				(((first_row_input & anding_mask)^(second_row_input & anding_mask))== 0)
					) 
					begin
						I_am_ready <=1;
					end
				else
					begin
						@(negedge clk);
						@(negedge clk);
						I_am_ready <=1;
					end
					

			end

		else 
			begin
				I_am_ready <=0;
			end

	end

				  

always @ (posedge clk)
	begin
		if(reset)
			begin
				ii <=0;	   
				iii<=0;	
			end
		else if(!reset) 
			begin
				if(ii < total/no_of_units && outsider5)
					begin
						package_by_package <= multipliers_output_vector;
						//@(posedge clk);
						ii <=ii+1;
					end
				else if(ii == total/no_of_units)
					begin
						package_by_package <= 0; 
					end
				
			end
	end							 

always @(posedge clk)
	begin

				if(iii <total/no_of_units -1)
					begin 
						if(final_adder_finish_dash) 
							begin
								iii <= iii+1;  
							end	 
							
					end
				else if(iii == total/no_of_units -1)
					begin 
						if(final_adder_finish_dash)
							begin	
								dot_product_output <= adder_output;
								finish<=1; 	  
							end
					end

	end




always @(posedge clk)
	begin  
		outsider1 <= outsider_read_now;	
		outsider2 <= outsider1;
		outsider3 <= outsider2 ;
		outsider4<=outsider3;
		outsider5<=outsider4; 

	end	


always @(posedge clk)
	begin
		if(outsider_read_now && flip)
			begin
				//demux_select <= 0;
				flip <= 1;
			end
		else if(/*!outsider_read_now && */(outsider1 || ~flip))
			begin
				if(flip)
					begin
					//	demux_select <= 1;
						first_row_four_elements_subset <= first_row_input[no_of_units*element_width-1-:(no_of_units/2)*element_width];
						second_row_four_elements_subset <= second_row_input[no_of_units*element_width-1-:(no_of_units/2)*element_width];
						flip <= 0;
					end
				else if(!flip)
					begin
					//	demux_select <=0;
						first_row_four_elements_subset <= first_row_input[(no_of_units/2)*element_width-1-:(no_of_units/2)*element_width];
						second_row_four_elements_subset <= second_row_input[(no_of_units/2)*element_width-1-:(no_of_units/2)*element_width];
						flip	<= 1;
					end
			end
	end	
	
	always @(posedge clk)
	begin
		if(outsider_read_now)
			begin
				demux_select <= 0;
				flip2 <= 1;
			end
		else if(!outsider_read_now && (outsider3 || ~flip2))
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