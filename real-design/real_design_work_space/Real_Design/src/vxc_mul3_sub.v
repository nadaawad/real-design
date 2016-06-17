module vXc_mul3_sub(clk,reset,first_row_plus_additional,constant,second_row_plus_additional,op,finish,counter2,counter5,result_mem_we,counter3,vXc_add_8_output); 

    parameter number_of_equations_per_cluster=16;
	parameter element_width_modified=34;
	parameter element_width=32;
	parameter no_of_units=8;
	parameter count=number_of_equations_per_cluster/no_of_units;
	parameter additional = no_of_units-(number_of_equations_per_cluster%no_of_units); 
	parameter total = number_of_equations_per_cluster+additional ;
	parameter number_of_clusters=1;
	integer counter=0;
	output reg[31:0] counter2;
	output reg[31:0] counter3;
	//integer counter4=0;
	output reg[31:0] counter5;
	integer i=0;
	
	input wire clk;
	input wire reset;		
	//input wire [element_width*number_of_equations_per_cluster-1:0]vector1;
	//input wire [element_width*number_of_equations_per_cluster-1:0]vector2;
	input wire [element_width-1:0] constant;
	input wire op;
	
	//output reg [element_width*number_of_equations_per_cluster-1:0]result;
			output reg finish;	
	output reg result_mem_we;//vector1_mem_we,vector2_mem_we;
	
	reg [no_of_units*element_width-1:0] first_row_input;
	reg [no_of_units*element_width-1:0] second_row_input;
	
	

	input wire [element_width*no_of_units-1:0] first_row_plus_additional;
		// pragma attribute first_row_plus_additional ram_block 1
	input wire [element_width*no_of_units-1:0] second_row_plus_additional;
		// pragma attribute second_row_plus_additional ram_block 1
	
		
	output wire [element_width*no_of_units-1:0]vXc_add_8_output;
	
	wire finish_out;
	
	
	vXc_add_8 #(.NI(no_of_units)) vXv(clk,reset,first_row_input,constant,second_row_input,op, vXc_add_8_output,finish_out );

	initial
		begin  
			counter3<=0; 
			counter2<=0;
			counter5<=0;   
			result_mem_we<=0;
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
					//		vector1_mem_we<=0;
			         //       vector2_mem_we<=0;
							//first_row_plus_additional[0] <= {vector1,{additional{32'b0}}};
							//second_row_plus_additional[0] <= {vector2,{additional{32'b0}}};
						end
						counter <= counter+1;
					end
				end
				
				always @ (posedge clk)
					begin
						if(reset)
							counter2<=0;
						else if(!reset)
							begin
								if(counter2 <total/no_of_units)
									begin
										if(counter2==0)
											first_row_input <= first_row_plus_additional;//[32*no_of_units*(total/no_of_units -(counter2-1)) -1 -:no_of_units*32];
										else
											begin
										@(posedge clk);
										@(posedge clk);
										first_row_input <= first_row_plus_additional;//[32*no_of_units*(total/no_of_units -(counter2-1)) -1 -:no_of_units*32];
										end
										counter2 <=counter2+1;
									end
								end
							end
							
							always @ (posedge clk)
								begin
									if(reset)
										counter5<=0;
									if(!reset)
										begin
											if(counter5 <total/no_of_units+2)
												begin  
													if(counter5==0)
														begin
															@(posedge clk);
															@(posedge clk);
															second_row_input <= second_row_plus_additional;//[32*no_of_units*(total/no_of_units -(counter5-1))-1 -:no_of_units*32];
														end
													else
														begin
															 @(posedge clk);
															@(posedge clk);
															second_row_input <= second_row_plus_additional;//[32*no_of_units*(total/no_of_units -(counter5-1))-1 -:no_of_units*32];
														end
														counter5 <=counter5+1;
													end
													end
												end
												
												always@(posedge clk)
													begin
														if(reset)
															begin
															counter3<=0;
															
															finish<=0;
															end
														if(finish_out)
															begin	
																
																if(counter3<total/no_of_units)
																	begin 
																		
																			
																		@(posedge clk);
																		result_mem_we<=1;
																		@(posedge clk);	
																		result_mem_we<=0;
																		/*if(op==1)
																			begin
																			@(posedge clk);
																		
																			@(posedge clk);	
																			end*/
																		
																		//result[element_width*(number_of_equations_per_cluster-(counter3-1)*no_of_units)-1-:no_of_units*element_width]<=vXc_add_8_output;
																		counter3<=counter3+1;
																	end
																else 
																	finish<=1;
																end
																end
										
									 
							endmodule