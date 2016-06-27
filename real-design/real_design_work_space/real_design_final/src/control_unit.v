
module control_unit (total,clk,reset,finish_alu,memories_pre_preprocess,memoryP_write_enable,memoryR_write_enable,memoryX_write_enable,memoryA_read_address,memoryP_read_address,memoryP_v2_read_address, memoryR_read_address,memoryX_read_address,memoryP_write_address,memoryR_write_address ,memoryX_write_address,halt,reset_vXv1,outsider_read_now,result_mem_we_4,memoryRprev_we,result_mem_we_5,result_mem_counter_5,read_again,start,read_again_2,result_mem_we_6,vXv1_finish,finish_all);
	
	parameter no_of_units = 8;
	parameter memory_read_address_width=32;	// m7tag yt3`yr
	parameter element_width = 32;
	
	integer counter=0;
	integer counter2=0;
	integer counter_vXv3=0;
	integer counter3=0;	
	integer counter4 = 0 ;	 
	integer counter5 = 0;	
	integer NumCyclesTillNow = 0;
	

	input wire clk,reset;
	input wire reset_vXv1;
	input wire outsider_read_now ;
	
	input wire result_mem_we_4;
	input wire read_again;	
	
	input wire result_mem_we_6;
	input wire read_again_2;	
	
	input wire finish_alu;
	
	input wire result_mem_we_5;
	input wire[31:0] result_mem_counter_5; 
	input wire start;
	input wire vXv1_finish;
	input wire finish_all;	
	input wire[31:0] total;
	input wire memories_pre_preprocess;
	
	
	
	output reg halt;
	reg increment_read_address_enable;
	reg [10:0]iteration_counter=0;
	
	
	
	output wire memoryX_write_enable; 
	assign memoryX_write_enable = result_mem_we_4;
	
	
	output wire memoryP_write_enable; 
	assign memoryP_write_enable = result_mem_we_6;
	
	
	
	//output reg memoryR_write_enable;
	output wire memoryR_write_enable;
	assign memoryR_write_enable=result_mem_we_5;
	
	output reg [memory_read_address_width-1:0]memoryA_read_address;
	output reg [memory_read_address_width-1:0]memoryP_read_address;	
	output reg [memory_read_address_width-1:0]memoryP_v2_read_address;
	output reg [memory_read_address_width-1:0]memoryX_read_address;
	output reg [memory_read_address_width-1:0]memoryR_read_address;
	output reg [memory_read_address_width-1:0]memoryX_write_address;
	output reg [memory_read_address_width-1:0]memoryP_write_address;
	//output reg [memory_read_address_width-1:0]memoryR_write_address;
	output wire [memory_read_address_width-1:0]memoryR_write_address;
	assign memoryR_write_address=result_mem_counter_5;
	output reg memoryRprev_we;
	
	
	reg finishvXv1_flag;	   // m7taga 23mlo reset fe kol iteration gdeda
	reg finish_start_flag;
	
	
	
	always @(posedge clk)
		begin
			if(reset)
			
				halt<=0;
		end
	
				
								
		
	/*	always@(posedge clk)
		begin 
			if(reset==1||finish_alu||memoryR_write_address>=(total/8))
				begin
					memoryR_write_address<=0;
					memoryR_write_enable<=0;
					
				end	
			
			
			else if(result_mem_we_5)
				begin 
					memoryR_write_enable<=1;
					memoryR_write_address<=result_mem_counter_5; 
				end	   

		
			else 
				begin
					memoryR_write_enable<=0;
					
					end
				
			end*/
			
	
			
			
			
			
			
			always@(posedge clk)
		
				begin
			
			
					if(reset==1||finish_alu)	    
				
						begin  
					
							finishvXv1_flag<=0;
					
							finish_start_flag<=0;
					
							memoryR_read_address<=0;
					
				
						end	
					
					
					else if(memoryR_read_address>= (total/8))
						
						begin
							memoryR_read_address<=0;
							finishvXv1_flag<=1;
							
							if(start)
								begin
									finish_start_flag<=1;	
								end
							end
			
			
			
			else if(!reset&&read_again_2&&!finish_alu)
				begin
					memoryR_read_address<= memoryR_read_address + 1'b1;
				end
				
				
			 else if(!reset_vXv1&&!finishvXv1_flag&&!finish_alu)
				begin	 
					memoryRprev_we<= 1; 
			
					@(posedge clk);
					
					if(!finishvXv1_flag)
						begin
		  			
							memoryR_read_address<= memoryR_read_address + 1'b1;
					  
							memoryRprev_we<= 1; 
					  end
				end	 
				
				
				 else if(start&&!finish_start_flag&&!finish_alu)
				begin	 
					 
					@(posedge clk);	
					if(!finish_start_flag)
		  		
						memoryR_read_address<= memoryR_read_address + 1'b1;
					   
				
					end	 
			
			   	end	 
		
		
			
		
		always@(posedge clk)
		begin
			if(reset||finish_alu)
				begin
					memoryP_v2_read_address <= 0;
			
					counter4 <=0 ;
					counter5 <=0 ;
				
				end
				
			else if (counter4 ==0 /* initialization*/ )	
			
				begin
			
					memoryP_v2_read_address <= 0;
			
					counter4 <=1 ;
			
				end	 
			
			else if((memoryP_v2_read_address >= (total/8)))	
			
				begin	
					memoryP_v2_read_address <= 0;
					counter5 <=1 ;
				end	
		
			else if(outsider_read_now && counter5==0)
				begin
					@(posedge clk);
					memoryP_v2_read_address	<= memoryP_v2_read_address + 1'b1;	
		   	
				end	
		
			else if	(read_again||read_again_2)
				begin 
				memoryP_v2_read_address	<= memoryP_v2_read_address + 1'b1;
				  end		
		
			end	 
		
		

	
	
		
			always@(posedge clk)
		
				begin	 
			
					NumCyclesTillNow = NumCyclesTillNow +1 ;
			
			
					if(reset==1||finish_alu)
				
						begin
					
							memoryA_read_address<=32'hffffffff;
					
				end
			
					else if(memories_pre_preprocess&&!halt)
				
						begin
					
							memoryA_read_address<=memoryA_read_address+1;
					
					end
			end	 
		
			
			
			always@(posedge clk)
	
				begin
			
					if(reset==1||finish_alu)//||memoryP_read_address>=number_of_clusters)
				
						begin
					
							memoryP_read_address<=0;
					
				
						end
			
					else if(increment_read_address_enable&&!halt)
				
						begin
					
							memoryP_read_address<=memoryP_read_address+1;  
					
					end
			end
			
			
			
			
			always@(posedge clk)
		
				begin
					if(reset==1||finish_alu||memoryP_write_address>=((total/8)))
						begin
							memoryP_write_address<=0;
					end
			
					else if(result_mem_we_6)
						begin
							memoryP_write_address<=memoryP_write_address+1;
						end
					end	
			
			
			
			
			
			
			
			always@(posedge clk)
		
				begin
			
					if(reset==1||finish_alu||memoryX_read_address>=((total/8)))
				
						begin
					
							memoryX_read_address<=0;
						end	
				
					else if(read_again)
		
						begin
	
							memoryX_read_address<= memoryX_read_address + 1'b1;	
		
						end
				
			end
			

			
		
			
			
			
		
			
			
			
			always @(posedge clk)
			begin
			
		if(reset==1||finish_alu||memoryX_write_address>=((total/8)) )
				begin
					memoryX_write_address<=0;
					
				end
			else if(result_mem_we_4)
				begin
					memoryX_write_address<=memoryX_write_address+1; 
					
				end
			end	  
			
			
			  											   
							
							
							
							
			
			
			
			
			
		
			
			
			
			always@(posedge clk)
			
				begin
				
					if(reset)
						begin
							counter3<=0;
							iteration_counter<=0;
							
						end
						
					else if(finish_all)
						
						begin
							
							iteration_counter<=iteration_counter+1;
							halt<=1;
						end
					else if(finish_alu&&!reset)
						
						begin
							counter3<=counter3+1;
							if(counter3==4)	 
								begin
									iteration_counter<=iteration_counter+1;

								end
							end
						else if(!finish_alu&&!reset)
							begin
								counter3<=0;
							end
						end
						
						
					
									
									
									
									
									
									
									
									
									
					
					//	
//					always@(posedge clk)
//						begin
//							$display("%d",memoryX_read_address);
//							
//							$display("%d",memoryX_write_address);
//								
//					end
//					

endmodule 
