//`define invalid 32'hFFFFFFFF   
`define invalid 24'hFFFFFF   // if we use 6 hex col now address 
`define smaller_large_zero 256'h0 
`define small_zero 32'h00000000

module P_Emap_8 (clk,input_data,read_preprocess,write_enable,write_address,col_nos,output_row,no_of_multiples,you_can_read,I_am_ready); 
	
	parameter no_of_elements_on_col_nos = 20 ;
	parameter no_of_elements_in_output = 8 ;
	parameter element_width = 32;   
	parameter no_of_units = 8;	
	parameter memory_height = 1000;
	parameter address_width= $clog2(memory_height)+1; 
	parameter col_nos_values_width = 32; 
	parameter multiples_memory_value_width = 32;
	
	
	input wire [multiples_memory_value_width-1:0] no_of_multiples; 	// this should come from the index matrix
	input wire I_am_ready;
	
	integer i = 1; 
	reg first_time_flag =1;	 
	output reg you_can_read=0;
	
input clk,read_preprocess,write_enable ; 
input wire [no_of_units*element_width-1:0]input_data; 
input wire [address_width - 1 : 0] write_address;
input [no_of_elements_on_col_nos*col_nos_values_width-1:0]col_nos;  // DONT MAKE THIS element_width
output wire [no_of_elements_in_output*element_width-1:0] output_row ;	



wire [col_nos_values_width-1 :0]   // // DONT MAKE THIS element_width
elem7_address,elem6_address,elem5_address,elem4_address,
elem3_address,elem2_address,elem1_address,elem0_address;


reg read_enable;


assign elem7_address = col_nos[((no_of_elements_on_col_nos-0-((i-1)*no_of_elements_in_output))) *col_nos_values_width-1-:col_nos_values_width];
assign elem6_address = col_nos[((no_of_elements_on_col_nos-1-((i-1)*no_of_elements_in_output))) *col_nos_values_width-1-:col_nos_values_width];
assign elem5_address = col_nos[((no_of_elements_on_col_nos-2-((i-1)*no_of_elements_in_output))) *col_nos_values_width-1-:col_nos_values_width];
assign elem4_address = col_nos[((no_of_elements_on_col_nos-3-((i-1)*no_of_elements_in_output))) *col_nos_values_width-1-:col_nos_values_width];
assign elem3_address = col_nos[((no_of_elements_on_col_nos-4-((i-1)*no_of_elements_in_output))) *col_nos_values_width-1-:col_nos_values_width];
assign elem2_address = col_nos[((no_of_elements_on_col_nos-5-((i-1)*no_of_elements_in_output))) *col_nos_values_width-1-:col_nos_values_width];
assign elem1_address = col_nos[((no_of_elements_on_col_nos-6-((i-1)*no_of_elements_in_output))) *col_nos_values_width-1-:col_nos_values_width];
assign elem0_address = col_nos[((no_of_elements_on_col_nos-7-((i-1)*no_of_elements_in_output))) *col_nos_values_width-1-:col_nos_values_width];

// (no_of_elements_on_col_nos-1-k-((i-1)*no_of_elements_in_output)) , k=0,1,2,3,.......,no_of_elements_in_output

reg [no_of_units * element_width - 1 : 0] mem [0 : memory_height];
// pragma attribute mem ram_block 1	

// d = devision , r = remainder
reg [col_nos_values_width-1:0 ]elem0_d;	reg [col_nos_values_width-1:0 ]elem0_r;  reg [col_nos_values_width-1:0 ]elem0_r_pipe; 
reg [col_nos_values_width-1:0 ]elem1_d;	reg [col_nos_values_width-1:0 ]elem1_r;  reg [col_nos_values_width-1:0 ]elem1_r_pipe; 
reg [col_nos_values_width-1:0 ]elem2_d;	reg [col_nos_values_width-1:0 ]elem2_r;  reg [col_nos_values_width-1:0 ]elem2_r_pipe;
reg [col_nos_values_width-1:0 ]elem3_d;	reg [col_nos_values_width-1:0 ]elem3_r;  reg [col_nos_values_width-1:0 ]elem3_r_pipe;
reg [col_nos_values_width-1:0 ]elem4_d;	reg [col_nos_values_width-1:0 ]elem4_r;  reg [col_nos_values_width-1:0 ]elem4_r_pipe;
reg [col_nos_values_width-1:0 ]elem5_d;	reg [col_nos_values_width-1:0 ]elem5_r;  reg [col_nos_values_width-1:0 ]elem5_r_pipe;
reg [col_nos_values_width-1:0 ]elem6_d;	reg [col_nos_values_width-1:0 ]elem6_r;  reg [col_nos_values_width-1:0 ]elem6_r_pipe;
reg [col_nos_values_width-1:0 ]elem7_d;	reg [col_nos_values_width-1:0 ]elem7_r;  reg [col_nos_values_width-1:0 ]elem7_r_pipe;



reg[no_of_units * element_width - 1 : 0] out7;
reg[no_of_units * element_width - 1 : 0] out6;
reg[no_of_units * element_width - 1 : 0] out5;
reg[no_of_units * element_width - 1 : 0] out4;
reg[no_of_units * element_width - 1 : 0] out3;
reg[no_of_units * element_width - 1 : 0] out2;
reg[no_of_units * element_width - 1 : 0] out1;
reg[no_of_units * element_width - 1 : 0] out0;   

assign output_row = {
out7[(no_of_units-elem7_r_pipe)*element_width-1-:element_width],
out6[(no_of_units-elem6_r_pipe)*element_width-1-:element_width],
out5[(no_of_units-elem5_r_pipe)*element_width-1-:element_width],
out4[(no_of_units-elem4_r_pipe)*element_width-1-:element_width],
out3[(no_of_units-elem3_r_pipe)*element_width-1-:element_width],
out2[(no_of_units-elem2_r_pipe)*element_width-1-:element_width],
out1[(no_of_units-elem1_r_pipe)*element_width-1-:element_width],
out0[(no_of_units-elem0_r_pipe)*element_width-1-:element_width]};

	
	initial 
		begin
			$readmemh("b.txt", mem);
		end	
		
	
		always @(posedge clk)
			begin 
				read_enable <= (read_preprocess || ~first_time_flag);
				
				
				//$display("%h",out15);  
				//$display("%h",out15[(no_of_units-elem15_r)*element_width-1-:element_width]); 
				//$display("%h",output_row);
			end	
		
	always @(posedge clk)
		begin 
			if(read_preprocess || ~first_time_flag)
				begin  
				if(elem0_address != `invalid)	
					begin elem0_d <= elem0_address/no_of_units; elem0_r <=elem0_address%no_of_units ;end
				else begin	elem0_d <= `invalid; elem0_r <=`small_zero; end	
				if(elem1_address != `invalid)	
					begin elem1_d <= elem1_address/no_of_units; elem1_r <=elem1_address%no_of_units ;end	   
				else begin	elem1_d <= `invalid; elem1_r <=`small_zero; end	
				if(elem2_address != `invalid)
					begin elem2_d <= elem2_address/no_of_units; elem2_r <=elem2_address%no_of_units ;end
				else begin	elem2_d <= `invalid; elem2_r <=`small_zero; end
				if(elem3_address != `invalid)
					begin elem3_d <= elem3_address/no_of_units; elem3_r <=elem3_address%no_of_units ;end
				else begin	elem3_d <= `invalid; elem3_r <=`small_zero; end
				if(elem4_address != `invalid)
					begin elem4_d <= elem4_address/no_of_units; elem4_r <=elem4_address%no_of_units ;end
				else begin	elem4_d <= `invalid; elem4_r <=`small_zero; end
				if(elem5_address != `invalid)
					begin elem5_d <= elem5_address/no_of_units; elem5_r <=elem5_address%no_of_units ;end	
				else begin	elem5_d <= `invalid; elem5_r <=`small_zero; end
				if(elem6_address != `invalid)
					begin elem6_d <= elem6_address/no_of_units; elem6_r <=elem6_address%no_of_units ;end
				else begin	elem6_d <= `invalid; elem6_r <=`small_zero; end
				if(elem7_address != `invalid)
					begin elem7_d <= elem7_address/no_of_units; elem7_r <=elem7_address%no_of_units ;end
				else begin	elem7_d <= `invalid; elem7_r <=`small_zero; end

		

				end	
			
		end	
	always @(posedge clk)
		begin  
		if(read_enable)	
			begin	 
				
				if(elem0_d != `invalid)	 
					out0<= mem[elem0_d];
				else out0<=`smaller_large_zero;	
				if(elem1_d != `invalid)	 
					out1<= mem[elem1_d]; 
				else out1<=`smaller_large_zero;	
				if(elem2_d != `invalid)	 
					out2<= mem[elem2_d];
				else out2<=`smaller_large_zero;	
				if(elem3_d != `invalid)	 
					out3<= mem[elem3_d];
				else out3<=`smaller_large_zero;	
				if(elem4_d != `invalid)	 
					out4<= mem[elem4_d];
				else out4<=`smaller_large_zero;	
				if(elem5_d != `invalid)	 
					out5<= mem[elem5_d];
				else out5<=`smaller_large_zero;	
				if(elem6_d != `invalid)	 
					out6<= mem[elem6_d];
				else out6<=`smaller_large_zero;	
				if(elem7_d != `invalid)	 
					out7<= mem[elem7_d]; 
				else out7<=`smaller_large_zero;	

			end 	
		end	   
		
	always @(posedge clk)
		begin
			if(read_enable)
				begin
					elem0_r_pipe <= elem0_r;
					elem1_r_pipe <= elem1_r;
					elem2_r_pipe <= elem2_r;
					elem3_r_pipe <= elem3_r;
					elem4_r_pipe <= elem4_r;
					elem5_r_pipe <= elem5_r;
					elem6_r_pipe <= elem6_r;
					elem7_r_pipe <= elem7_r;
					
				end	
			
		end	  
		
		always @(posedge clk)
			begin
				if(read_preprocess || ~first_time_flag)	
					begin 
						you_can_read<=1;
						@(posedge clk);	
						you_can_read<=0;

						if(i<no_of_multiples)
							begin
								i<=i+1;
								first_time_flag<=0;
							end	 
						else if(i>= no_of_multiples)
							begin 
								i<=1; 
								first_time_flag<=1;
							end	   
					end	 
				else you_can_read<=0;	
			end	  
			
			
				
			always @(posedge clk) 
		begin
			if( write_enable == 1'b1 ) 
				begin
					mem[write_address] <= input_data; 
				end
			end

endmodule 	 
	 