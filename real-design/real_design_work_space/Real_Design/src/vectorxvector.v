module vectorXvector(clk,reset,first_row_plus_additional,second_row_plus_additional,result,finish); 

    parameter number_of_equations_per_cluster=10;
	parameter element_width_modified=34;
	parameter element_width=32;
	parameter no_of_units=8;
	parameter count=number_of_equations_per_cluster/no_of_units;
	parameter additional = no_of_units-(number_of_equations_per_cluster%no_of_units); 
	parameter total = number_of_equations_per_cluster+additional ;
	parameter number_of_clusters=1;
	integer counter ;
	integer counter2;
	
	integer i=0;
	
	input wire clk;
	input wire reset;		
	//input wire [element_width*number_of_equations_per_cluster-1:0]vector1;
//	input wire [element_width*number_of_equations_per_cluster-1:0]vector2;
	output wire [element_width-1:0]result;
	output wire finish;
	
	
    
	
	input wire [element_width*(no_of_units)-1:0] first_row_plus_additional;
	input wire [element_width*(no_of_units)-1:0] second_row_plus_additional;	   
	
	
	
	reg [no_of_units*element_width-1:0] first_row_input;
	reg [no_of_units*element_width-1:0] second_row_input;
	
	//
//	reg [element_width*(total)-1:0] first_row_plus_additional[0:0];
//	// pragma attribute first_row_plus_additional ram_block 1
//	reg [element_width*(total)-1:0] second_row_plus_additional[0:0];
//		// pragma attribute second_row_plus_additional ram_block 1
	
		
		reg [element_width-1:0]dot_product_output;
	
	
	eight_Dot_Product_Multiply #(.NOE(number_of_equations_per_cluster))
	vXv(clk,reset,first_row_input,second_row_input, result,finish );

//	Sixteen_Dot_Product_Multiply #(.NOE(number_of_equations_per_cluster))
	// vXv(clk,reset,first_row_input,second_row_input, result,finish );

//onezerotwofour_Dot_Product_Multiply #(.NOE(number_of_equations_per_cluster))
	//vXv(clk,reset,first_row_input,second_row_input, result,finish );
	always@(posedge clk)
		begin
			if(reset)
				begin
					counter<=0;	
				    
                   
				end
			else
				begin
					if (counter==0)
						begin
							
							counter <= counter+1;
						end
					
						

				end
		end

	always @ (posedge clk)
					begin
					if(reset)
						begin
						counter2<=0;
						end
						else 
							begin
								@(posedge clk);
								if(counter2 <total/no_of_units+2)
									begin  
										first_row_input <= first_row_plus_additional;
										second_row_input <= second_row_plus_additional;
										counter2 <=counter2+1;
									end
								end
							end

								
							endmodule