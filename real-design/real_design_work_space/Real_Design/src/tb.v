`timescale 1ns / 1ps

module tb();

	reg clk;
	reg reset;
	wire halt;
	wire finish;
	reg reset_vXv1;
	reg reset_mXv1;
	
	integer counter1= 0;
	integer counter2=0;	
	integer counter3=0;
	
	top_module uut (clk,reset,finish,reset_vXv1,reset_mXv1,halt);
	
	reg[31:0] AA = 32'h42c80000;
	reg[31:0] BB = 32'h42c80000;
	wire [31:0]CC ;	  
	wire [31:0]DD ;	
	wire [31:0]EE ;
	
	

	initial
		begin
			clk<=0;
			reset <= 1;
			reset_vXv1<=1;
			reset_mXv1<=0;
		end
		
		
		
		always @(posedge clk)
			begin
				
				if(!finish)	 
					
					begin
						 counter3<=0;
						if(counter1==0)
							begin
								counter1= counter1+1;
							end
						else if(counter1 ==1)
							begin
								reset <=0;
								counter1 = counter1 +1;
							end
						else if(counter1==2)
							begin
								reset_mXv1<=1;
								reset_vXv1<=0; 
								counter1 = counter1+1;
							end
						end
						
						
						else
							begin
								reset_vXv1<=1;
								reset_mXv1<=0;
								counter3<=counter3+1;
								if(counter3==4)
									
									begin
							    
										reset_vXv1<=0;
								
										reset_mXv1<=1;
								end
								
							end 
						end
						
						always 
							begin
								#10 clk=~clk;
							end
							
							always @(posedge clk)
								begin
									if(halt == 1)
										begin
											$display (halt);
											$finish();
										end
									end
								endmodule
