module vectorXvector_mXv(clk,reset,first_row_plus_additional,vector2,result,finish,AP_total_mem_we,counter2); 

    parameter number_of_equations_per_cluster=16;
	parameter element_width_modified=34;
	parameter element_width=32;
	parameter no_of_units=8;
	parameter count=number_of_equations_per_cluster/no_of_units;
	parameter additional = no_of_units-(number_of_equations_per_cluster%no_of_units); 
	parameter total = number_of_equations_per_cluster+additional ;
	parameter number_of_clusters=1;
	integer counter = 0;
	integer i=0;
	integer counter3=0;
	

	input wire clk;
	input wire reset;		
	input wire [element_width*no_of_units-1:0]vector2;
	input wire [element_width*no_of_units:0] first_row_plus_additional;
	//input wire [element_width*total-1:0] AP_total ;
	
	
	
	//output wire [element_width*number_of_equations_per_cluster-1:0] AP;
	output wire [element_width-1:0]result;
	output wire finish;		 
	
	output reg AP_total_mem_we;	 
	output reg[31:0] counter2;
	
	reg [no_of_units*element_width-1:0] first_row_input;
	reg [no_of_units*element_width-1:0] second_row_input;
	reg [element_width-1:0]dot_product_output;
	
	
	
	
	//assign AP=AP_total[element_width*total-1:element_width*total-element_width*number_of_equations_per_cluster];
	

	eight_Dot_Product_Multiply #(.NOE(number_of_equations_per_cluster))
	vXv(clk,reset,first_row_input,second_row_input, result,finish );
		
//fiveonetwo_Dot_Product_Multiply#(.NOE(number_of_equations_per_cluster))
//vXv(clk,reset,first_row_input,second_row_input, result,finish );

//Sixteen_Dot_Product_Multiply #(.NOE(number_of_equations_per_cluster))
//vXv(clk,reset,first_row_input,second_row_input, result,finish );

//onezerotwofour_Dot_Product_Multiply #(.NOE(number_of_equations_per_cluster))
	//vXv(clk,reset,first_row_input,second_row_input, result,finish );
	
	initial
		begin
		
		  counter2<=0;
		end
		



	always @(posedge clk)
		begin
			if(reset)
				begin
				counter<=0;
				
				end
			else if(!reset)
				begin
					if (counter==0)
						begin 
							
							//first_row_plus_additional[0] <= {vector1,{additional{32'b0}}};
						end
						counter <= counter+1;
					end
				end
				
				always @ (posedge clk)
					begin 
						if(reset)
							begin
							counter2<=0;
							AP_total_mem_we<=0;
							end
						else if(!reset)
							begin  

								if(counter2 <total/no_of_units)
									begin 
										AP_total_mem_we<=1;
										@(posedge clk);
										AP_total_mem_we<=0;
										first_row_input <= first_row_plus_additional;
										second_row_input <= vector2; 
										counter2 <=counter2+1;
									end	 
								
								end
							end
							
							 
							endmodule