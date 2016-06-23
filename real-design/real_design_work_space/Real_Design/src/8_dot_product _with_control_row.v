
module eight_Dot_Product_Multiply_with_control_row(clk,reset ,first_row_input,second_row_input, dot_product_output,finish,outsider_read_now,no_of_multiples,prepare_my_new_input,fake_prepare0,I_am_ready);

parameter NOE = 10;
parameter NI = 8;
integer repetition_times = (NI==8)?6:(NI==16)?6:0;
parameter additional = NI-(NOE%NI); 
parameter total = NOE+additional ;	 

reg initialization_counter = 1;
reg initialization_counter2 = 1;  


integer ii=0; 
integer iii=0; 
integer iiii=0;
integer iiiii=0;
output reg finish ;	 
output reg prepare_my_new_input=0; 

output reg fake_prepare0=0;
reg fake_prepare1=0;
reg fake_prepare2=0; 

reg fake_reset;

input wire [31:0] no_of_multiples;
reg second_number_of_multiples = 10000;
reg [31:0] delayed_no_of_multiples=10000;

input wire outsider_read_now;
input wire reset ;
input wire[32*NI-1:0] first_row_input;
input wire[32*NI-1:0] second_row_input;

reg [32*NI-1:0] anding_mask = {NI{32'h80000000}};  
output reg I_am_ready;


reg save = 0;
reg adder_tree_start=0;
input clk ;
reg[32*NI-1:0] package_by_package;
wire [32*NI-1:0] multipliers_output_vector;
wire [31:0] adder_output;
output reg [31:0] dot_product_output;


wire[32*(NI/2)-1:0] demux_four_inputs;
reg demux_select;
reg flip=1; 
reg flip2=1;

//four_to_eight_demux demux_1(demux_four_inputs,demux_select,multipliers_output_vector);
N_to_2N_demux  #(.NI(NI))demux(demux_four_inputs,demux_select,multipliers_output_vector);
reg [32*(NI/2)-1:0] first_row_four_elements_subset;
reg [32*(NI/2)-1:0] second_row_four_elements_subset;	  

reg outsider1=0;
reg outsider2=0;
reg outsider3=0;
reg outsider4=0;  
reg outsider5=0;   


reg reset_pip1,reset_pip2,reset_pip3;

reg[3:0] fifo_read_address =1;
reg[3:0] fifo_write_address=0 ;
reg fifo_write_enable=0; 

wire[31:0] fifo_output_data;

wire ExE_finish;
wire final_adder_finish_dash;

genvar j ;
generate
for(j=0;j<NI/2;j=j+1) begin : instantiate_Multiplier

multiply m (first_row_four_elements_subset[32*(NI/2-j)-1-:32], second_row_four_elements_subset[32*(NI/2-j)-1-:32], clk, 1, demux_four_inputs[32*(NI/2-j)-1-:32]);
end
endgenerate


Eight_Organizer_with_control_row #(.NI(NI)) E_O (clk,package_by_package,adder_tree_start , adder_output,outsider4,final_adder_finish_dash,ExE_finish);

multiples_fifo M_F (clk,fifo_write_enable,fifo_read_address,fifo_write_address,no_of_multiples,fifo_output_data);




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

				  


always@(posedge clk)
	begin 
		if(reset)
			begin
				prepare_my_new_input<=0;	
			end	   
		
		else if(iiiii <(no_of_multiples-1) && outsider3) 
			begin
				iiiii <=iiiii+1;
			end
		else if(iiiii==(no_of_multiples-1) && outsider3)
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
 
		if(outsider_read_now && ! initialization_counter)
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
		if(outsider_read_now && initialization_counter)	
			begin 
				delayed_no_of_multiples <= no_of_multiples;
				initialization_counter<=0;
			end	 
		else if(fake_reset)
			begin
				delayed_no_of_multiples <=fifo_output_data;
				fifo_read_address <= (fifo_read_address +1)%8;


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
			if(ii < delayed_no_of_multiples && outsider5)
				begin
					package_by_package <= multipliers_output_vector;
					//@(posedge clk);
					ii <=ii+1;
				end
			else if(ii == delayed_no_of_multiples)
				begin
					package_by_package <= 0; 
				end
			
		end
end							 

always @(posedge clk)
	begin

				if(iii <delayed_no_of_multiples -1)
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
		if(outsider_read_now)
			begin
				//demux_select <= 0;
				flip <= 1;
			end
		else if(!outsider_read_now && (outsider1 || ~flip))
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