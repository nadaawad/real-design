module vectorXvector_mXv_with_control(total,clk,reset,first_row_plus_additional,vector2,result,finish,AP_total_mem_we,counter2,outsider_read_now); 

   
	parameter element_width=32;
	parameter no_of_units=8;
	
	integer counter = 0;
	integer i=0;
	integer counter3=0;
	
	
	input wire outsider_read_now;
	input wire clk;
	input wire reset;		
	input wire [element_width*no_of_units-1:0]vector2;
	input wire [element_width*no_of_units-1:0] first_row_plus_additional;
	input wire [31:0]total;
	
	
	
	//output wire [element_width*number_of_equations_per_cluster-1:0] AP;
	output wire [element_width-1:0]result;
	output wire finish;		 
	
	output reg AP_total_mem_we=0;	 
	output reg[31:0] counter2;
	
	reg [no_of_units*element_width-1:0] first_row_input;
	reg [no_of_units*element_width-1:0] second_row_input;
	
	
	

	eight_Dot_Product_Multiply_with_control #(.no_of_units(no_of_units),.element_width(element_width))
	vXv(clk,reset ,first_row_input,second_row_input, dot_product_output,finish,outsider_read_now,total,I_am_ready);
		

	
	initial
		begin
		
		  counter2<=0;
		end
		




				
				always @ (posedge clk)
					begin 
						if(reset)
							begin
							counter2<=0;
						
							end
						else if(!reset)
							begin  
								if(counter2 <total/no_of_units+2 && outsider_read_now)
									begin 
										
										first_row_input <= first_row_plus_additional;
										second_row_input <= vector2; 
										
										AP_total_mem_we<=1;
										@(posedge clk);
										AP_total_mem_we<=0;
										counter2 <=counter2+1;
									end	 
								
								end
							end
							
							 
							endmodule