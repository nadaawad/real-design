`define invalid 32'hFFFFFFFF  
`define large_zero 512'h00000000 
`define small_zero 32'h00000000

module P_Emap_16 (clk,read_preprocess,write_enable,col_nos,output_row,no_of_multiples); 
	
	parameter no_of_elements_on_col_nos = 32 ;
	parameter no_of_elements_in_output = 16 ;
	parameter element_width = 32;   
	parameter no_of_units = 8;	
	input wire [31:0] no_of_multiples; 	// this should come from the index matrix

	
	integer i = 1; 
	reg first_time_flag =1;
	
input clk,read_preprocess,write_enable ;
input [no_of_elements_on_col_nos*element_width-1:0]col_nos;
output wire [no_of_elements_in_output*element_width-1:0] output_row ;	



wire [element_width-1 :0]
elem15_address,elem14_address,elem13_address,elem12_address,
elem11_address,elem10_address,elem9_address,elem8_address,
elem7_address,elem6_address,elem5_address,elem4_address,
elem3_address,elem2_address,elem1_address,elem0_address;


reg read_enable;

assign elem15_address = col_nos[((no_of_elements_on_col_nos-0-((i-1)*no_of_elements_in_output)))*element_width-1-:element_width];
assign elem14_address = col_nos[((no_of_elements_on_col_nos-1-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];	
assign elem13_address = col_nos[((no_of_elements_on_col_nos-2-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem12_address = col_nos[((no_of_elements_on_col_nos-3-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem11_address = col_nos[((no_of_elements_on_col_nos-4-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem10_address = col_nos[((no_of_elements_on_col_nos-5-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem9_address = col_nos[((no_of_elements_on_col_nos-6-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem8_address = col_nos[((no_of_elements_on_col_nos-7-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem7_address = col_nos[((no_of_elements_on_col_nos-8-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem6_address = col_nos[((no_of_elements_on_col_nos-9-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem5_address = col_nos[((no_of_elements_on_col_nos-10-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem4_address = col_nos[((no_of_elements_on_col_nos-11-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem3_address = col_nos[((no_of_elements_on_col_nos-12-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem2_address = col_nos[((no_of_elements_on_col_nos-13-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem1_address = col_nos[((no_of_elements_on_col_nos-14-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];
assign elem0_address = col_nos[((no_of_elements_on_col_nos-15-((i-1)*no_of_elements_in_output))) *element_width-1-:element_width];

// (no_of_elements_on_col_nos-1-k-((i-1)*no_of_elements_in_output)) , k=0,1,2,3,.......,no_of_elements_in_output

reg [no_of_units * element_width - 1 : 0] mem [0 : 100000];
// pragma attribute mem ram_block 1	

// d = devision , r = remainder
reg [element_width-1:0 ]elem0_d;	reg [element_width-1:0 ]elem0_r;  reg [element_width-1:0 ]elem0_r_pipe; 
reg [element_width-1:0 ]elem1_d;	reg [element_width-1:0 ]elem1_r;  reg [element_width-1:0 ]elem1_r_pipe; 
reg [element_width-1:0 ]elem2_d;	reg [element_width-1:0 ]elem2_r;  reg [element_width-1:0 ]elem2_r_pipe;
reg [element_width-1:0 ]elem3_d;	reg [element_width-1:0 ]elem3_r;  reg [element_width-1:0 ]elem3_r_pipe;
reg [element_width-1:0 ]elem4_d;	reg [element_width-1:0 ]elem4_r;  reg [element_width-1:0 ]elem4_r_pipe;
reg [element_width-1:0 ]elem5_d;	reg [element_width-1:0 ]elem5_r;  reg [element_width-1:0 ]elem5_r_pipe;
reg [element_width-1:0 ]elem6_d;	reg [element_width-1:0 ]elem6_r;  reg [element_width-1:0 ]elem6_r_pipe;
reg [element_width-1:0 ]elem7_d;	reg [element_width-1:0 ]elem7_r;  reg [element_width-1:0 ]elem7_r_pipe;
reg [element_width-1:0 ]elem8_d;	reg [element_width-1:0 ]elem8_r;  reg [element_width-1:0 ]elem8_r_pipe;
reg [element_width-1:0 ]elem9_d;	reg [element_width-1:0 ]elem9_r;  reg [element_width-1:0 ]elem9_r_pipe;
reg [element_width-1:0 ]elem10_d;	reg [element_width-1:0 ]elem10_r; reg [element_width-1:0 ]elem10_r_pipe;
reg [element_width-1:0 ]elem11_d;	reg [element_width-1:0 ]elem11_r; reg [element_width-1:0 ]elem11_r_pipe;
reg [element_width-1:0 ]elem12_d;	reg [element_width-1:0 ]elem12_r; reg [element_width-1:0 ]elem12_r_pipe;
reg [element_width-1:0 ]elem13_d;	reg [element_width-1:0 ]elem13_r; reg [element_width-1:0 ]elem13_r_pipe;
reg [element_width-1:0 ]elem14_d;	reg [element_width-1:0 ]elem14_r; reg [element_width-1:0 ]elem14_r_pipe;
reg [element_width-1:0 ]elem15_d;	reg [element_width-1:0 ]elem15_r; reg [element_width-1:0 ]elem15_r_pipe;


reg[no_of_units * element_width - 1 : 0] out15;
reg[no_of_units * element_width - 1 : 0] out14;
reg[no_of_units * element_width - 1 : 0] out13;
reg[no_of_units * element_width - 1 : 0] out12;
reg[no_of_units * element_width - 1 : 0] out11;
reg[no_of_units * element_width - 1 : 0] out10;
reg[no_of_units * element_width - 1 : 0] out9;
reg[no_of_units * element_width - 1 : 0] out8;
reg[no_of_units * element_width - 1 : 0] out7;
reg[no_of_units * element_width - 1 : 0] out6;
reg[no_of_units * element_width - 1 : 0] out5;
reg[no_of_units * element_width - 1 : 0] out4;
reg[no_of_units * element_width - 1 : 0] out3;
reg[no_of_units * element_width - 1 : 0] out2;
reg[no_of_units * element_width - 1 : 0] out1;
reg[no_of_units * element_width - 1 : 0] out0;   

assign output_row = {out15[(no_of_units-elem15_r_pipe)*element_width-1-:element_width],
out14[(no_of_units-elem14_r_pipe)*element_width-1-:element_width],
out13[(no_of_units-elem13_r_pipe)*element_width-1-:element_width],
out12[(no_of_units-elem12_r_pipe)*element_width-1-:element_width],
out11[(no_of_units-elem11_r_pipe)*element_width-1-:element_width],
out10[(no_of_units-elem10_r_pipe)*element_width-1-:element_width],
out9[(no_of_units-elem9_r_pipe)*element_width-1-:element_width],
out8[(no_of_units-elem8_r_pipe)*element_width-1-:element_width],
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
				$display("%h",output_row);
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
				if(elem8_address != `invalid)
					begin elem8_d <= elem8_address/no_of_units; elem8_r <=elem8_address%no_of_units ;end	
				else begin	elem8_d <= `invalid; elem8_r <=`small_zero; end
				if(elem9_address != `invalid)
					begin elem9_d <= elem9_address/no_of_units; elem9_r <=elem9_address%no_of_units ;end
				else begin	elem9_d <= `invalid; elem9_r <=`small_zero; end
				if(elem10_address != `invalid)
					begin elem10_d <= elem10_address/no_of_units; elem10_r <=elem10_address%no_of_units ;end
				else begin	elem10_d <= `invalid; elem10_r <=`small_zero; end
				if(elem11_address != `invalid)
					begin elem11_d <= elem11_address/no_of_units; elem11_r <=elem11_address%no_of_units ;end
				else begin	elem11_d <= `invalid; elem11_r <=`small_zero; end
				if(elem12_address != `invalid)
					begin elem12_d <= elem12_address/no_of_units; elem12_r <=elem12_address%no_of_units ;end
				else begin	elem12_d <= `invalid; elem12_r <=`small_zero; end
				if(elem13_address != `invalid)
					begin elem13_d <= elem13_address/no_of_units; elem13_r <=elem13_address%no_of_units ;end
				else begin	elem13_d <= `invalid; elem13_r <=`small_zero; end
				if(elem14_address != `invalid)
					begin elem14_d <= elem14_address/no_of_units; elem14_r <=elem14_address%no_of_units ;end
				else begin	elem14_d <= `invalid; elem14_r <=`small_zero; end
				if(elem15_address != `invalid)
					begin elem15_d <= elem15_address/no_of_units; elem15_r <=elem15_address%no_of_units ;end
				else begin	elem15_d <= `invalid; elem15_r <=`small_zero; end

				end	
			
		end	
	always @(posedge clk)
		begin  
		if(read_enable)	
			begin	 
				
				if(elem0_d != `invalid)	 
					out0<= mem[elem0_d];
				else out0<=`large_zero;	
				if(elem1_d != `invalid)	 
					out1<= mem[elem1_d]; 
				else out1<=`large_zero;	
				if(elem2_d != `invalid)	 
					out2<= mem[elem2_d];
				else out2<=`large_zero;	
				if(elem3_d != `invalid)	 
					out3<= mem[elem3_d];
				else out3<=`large_zero;	
				if(elem4_d != `invalid)	 
					out4<= mem[elem4_d];
				else out4<=`large_zero;	
				if(elem5_d != `invalid)	 
					out5<= mem[elem5_d];
				else out5<=`large_zero;	
				if(elem6_d != `invalid)	 
					out6<= mem[elem6_d];
				else out6<=`large_zero;	
				if(elem7_d != `invalid)	 
					out7<= mem[elem7_d]; 
				else out7<=`large_zero;	
				if(elem8_d != `invalid)	 
					out8<= mem[elem8_d];
				else out8<=`large_zero;	
				if(elem9_d != `invalid)	 
					out9<= mem[elem9_d];
				else out9<=`large_zero;	
				if(elem10_d != `invalid)	 
					out10<= mem[elem10_d];
				else out10<=`large_zero;	
				if(elem11_d != `invalid)	 
					out11<= mem[elem11_d];
				else out11<=`large_zero;	
				if(elem12_d != `invalid)	 
					out12<= mem[elem12_d]; 
				else out12<=`large_zero;	
				if(elem13_d != `invalid)	 
					out13<= mem[elem13_d]; 
				else out13<=`large_zero;	
				if(elem14_d != `invalid)	 
					out14<= mem[elem14_d];
				else out14<=`large_zero;	
				if(elem15_d != `invalid)	 
					out15<= mem[elem15_d];
				else out15<=`large_zero;	
					
 
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
					elem8_r_pipe <= elem8_r;
					elem9_r_pipe <= elem9_r;
					elem10_r_pipe <= elem10_r;
					elem11_r_pipe <= elem11_r;
					elem12_r_pipe <= elem12_r;
					elem13_r_pipe <= elem13_r;
					elem14_r_pipe <= elem14_r;
					elem15_r_pipe <= elem15_r;
					
				end	
			
		end	  
		
		always @(posedge clk)
			begin
				if(read_preprocess || ~first_time_flag)	
					begin
						if(i<no_of_multiples)
							begin
								i<=i+1;
								first_time_flag<=0;
							end	 
						else
							begin 
								i<=1; 
								first_time_flag<=1;
							end	   
					end	
			end	

endmodule 	 
	 