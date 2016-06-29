
module eight_Dot_Product_Multiply_with_control_row(clk,main_reset,reset ,first_row_input,second_row_input, dot_product_output,finish,outsider_read_now,no_of_multiples,prepare_my_new_input,fake_prepare0,I_am_ready);


parameter no_of_units = 8;
parameter element_width = 32;
parameter multiples_memory_value_width = 3;

reg initialization_counter = 1;
reg pos_neg_indicator = 0;  



integer ii=0; 
integer iii=0; 
integer iiii=0;
integer iiiii=0;
output reg finish=0 ;	 
output reg prepare_my_new_input=0; 

output reg fake_prepare0=0;
reg fake_prepare1=0;
reg fake_prepare2=0; 

reg fake_reset=0;

input wire main_reset;
input wire [multiples_memory_value_width-1:0] no_of_multiples;
reg [multiples_memory_value_width-1:0] delayed_no_of_multiples=10000;

input wire outsider_read_now;
input wire reset ;
input wire[element_width*no_of_units-1:0] first_row_input;
input wire[element_width*no_of_units-1:0] second_row_input;

reg [element_width*no_of_units-1:0] anding_mask = {no_of_units{32'h80000000}};  // IF 64 CHANGE MANULAY
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


reg reset_pip1,reset_pip2,reset_pip3;

reg[3:0] fifo_read_address =0;
reg[3:0] fifo_write_address=0 ;
reg fifo_write_enable=0; 

wire[multiples_memory_value_width-1:0] fifo_output_data;

wire ExE_finish;
wire final_adder_finish_dash;

genvar j ;
generate
for(j=0;j<no_of_units/2;j=j+1) begin : instantiate_Multiplier

multiply m (first_row_four_elements_subset[element_width*(no_of_units/2-j)-1-:element_width], second_row_four_elements_subset[element_width*(no_of_units/2-j)-1-:element_width], clk, 1, demux_four_inputs[element_width*(no_of_units/2-j)-1-:element_width]);
end
endgenerate


Eight_Organizer_with_control_row #(.no_of_units(no_of_units)) E_O (clk,package_by_package,adder_tree_start , adder_output,outsider4,final_adder_finish_dash,ExE_finish,main_reset);

multiples_fifo #(.multiples_memory_value_width(multiples_memory_value_width)) M_F (clk,fifo_write_enable,fifo_read_address,fifo_write_address,no_of_multiples,fifo_output_data);


initial 
begin
initialization_counter <= 1;
end

always@(posedge clk)
	begin 
		if(reset)
			begin
				prepare_my_new_input<=0;	
			end	   
		
		else if(iiiii <(no_of_multiples-1) && outsider2) 
			begin
				iiiii <=iiiii+1;
			end
		else if(iiiii==(no_of_multiples-1) && outsider2)
		    begin 
				prepare_my_new_input<=1;

			end	
		 else if(iiiii >=(no_of_multiples-1) && prepare_my_new_input)
			begin
				iiiii <=0; 	   // THIS HAS TO BE CHANGED 
				prepare_my_new_input<=0;
			end	
		  
			
		 if(iiii <(delayed_no_of_multiples-1) && ExE_finish) 
			begin
				iiii <=iiii+1;
			end	 		
		 else if(iiii==(delayed_no_of_multiples-1) && ExE_finish)
		    begin 
				fake_prepare0<=1;
			end	 
 
		 else if(iiii >=(delayed_no_of_multiples-1) && fake_prepare0)
			begin
				iiii <=0; 	   // THIS HAS TO BE CHANGED 
				fake_prepare0<=0;
			end	
	end	



always @(posedge clk)
	begin 
		 if(reset )
			begin
				fifo_write_address <= (fifo_write_address +1)%10 ;
				fifo_write_enable<=1;	

			end
		else
			begin
				fifo_write_enable<=0;
			end
	end	


always @(posedge clk)
	begin
		if(main_reset)
		 begin initialization_counter <=1;end
		else if(reset_pip2 && initialization_counter)
		begin
				delayed_no_of_multiples <=fifo_output_data;
				fifo_read_address <= (fifo_read_address +1)%10;
				initialization_counter<=0;
		end
		
		else if(fake_reset)
			begin
				delayed_no_of_multiples <=fifo_output_data;
				fifo_read_address <= (fifo_read_address +1)%10;
			end
	end	

always @ (posedge clk)
begin

 if(!reset) 
		begin
			if(ii < delayed_no_of_multiples && outsider5)
				begin
					package_by_package <= multipliers_output_vector;
					//@(posedge clk);
					ii <=ii+1;
				end
			else if(ii == delayed_no_of_multiples)
				begin
					package_by_package <= 0; 
					ii<=0;
				end
			
		end
end							 

always @(posedge clk)
	begin

				if(main_reset)
				begin finish <=0 ; end 

				else if(iii <delayed_no_of_multiples -1)
					begin 
						if(final_adder_finish_dash) 
							begin
								iii <= iii+1;  
							end	 
							
					end
				else if(iii == delayed_no_of_multiples -1)
					begin 
						if(final_adder_finish_dash)
							begin	
								dot_product_output <= adder_output;
								finish<=1; 	 
								iii<=0;
							end
					end

	end




always @(posedge clk)
begin
	if(fake_prepare2) 
		begin 
			adder_tree_start <= 0;
		end
	else if(!fake_prepare2)
	begin
	//if(outsider_read_now) begin adder_tree_start <=1; end	
		adder_tree_start <=1;
	end
end

always @(posedge clk)
	begin
	// $display("DOT PRODUCT OUTPU : %h",dot_product_output);	
	end	

always @(posedge clk)
	begin  
		outsider1 <= outsider_read_now;	
		outsider2 <= outsider1;
		outsider3 <= outsider2 ;
		outsider4<=outsider3;
		outsider5<=outsider4; 

		
	    fake_prepare1<=fake_prepare0;
		fake_prepare2<=fake_prepare1; 
		fake_reset<=fake_prepare2;

		reset_pip1 <=reset;
		reset_pip2 <=reset_pip1;	
		reset_pip3 <=reset_pip2;	
	
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