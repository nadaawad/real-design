
`timescale 1 ns / 1 ps

							 
	
`define zero_filling 32'd0
module matrix_by_vector_v3_with_control (clk,reset,start,A_rows,vector_rows,out,finish,outsider_read_now,no_of_multiples,total_with_additional_A,memories_pre_preprocess,you_can_read,I_am_ready);  

	// NOTE THAT : number of multiples determines how many cycles should the row by vector work on a certain row before it finishes
	// While total_with_additional_A / no_of_row_by_vector_modules determines the number of times we will give out row_by_vector modules
	// new rows to operate on .
	
parameter element_width =32;   
parameter no_of_row_by_vector_modules =4; 
parameter NI = 8;

         

//}} End of automatically maintained section

// -- Enter your statements here -- //		
input clk,reset,start;	 
input wire [no_of_row_by_vector_modules*element_width-1:0] no_of_multiples;

wire clk,reset,start;
input A_rows;
wire[no_of_row_by_vector_modules*element_width*NI-1:0] A_rows;  
input vector_rows;
wire[no_of_row_by_vector_modules*NI*element_width-1:0]vector_rows;	  
input[31:0] total_with_additional_A;

input wire[no_of_row_by_vector_modules-1:0] you_can_read;

wire[element_width*(no_of_row_by_vector_modules)-1:0] result;
reg sel;

output out;
wire[NI*element_width-1:0] out;

output reg finish;

output reg memories_pre_preprocess=0;

integer i=0;
integer counter=0;
integer counter2=0;	




wire [no_of_row_by_vector_modules-1:0] give_us_all;
reg [no_of_row_by_vector_modules-1:0] give_us_all_dash=0;

reg first_initialization_counter=0;

reg first_pipeline=1;	
reg second_pipeline=1;


wire[no_of_row_by_vector_modules-1:0] decoder_read_now;	 
reg [no_of_row_by_vector_modules-1:0] decoder_read_now_dash=0;
reg [no_of_row_by_vector_modules-1:0] start_row_by_vector;
output wire [no_of_row_by_vector_modules-1:0] I_am_ready;	  

reg matrix_by_vector_finished ;	 
output wire outsider_read_now;		


// DECODER 

decoder_with_control #(.no_of_row_by_vector_modules(no_of_row_by_vector_modules)) d (clk,result ,out , (&decoder_read_now_dash),outsider_read_now);

//



genvar j;
generate
for(j=0;j<no_of_row_by_vector_modules;j=j+1) begin:instantiate_ROW_BY_VECTOR
	
row_by_vector_with_control #(.NI(NI),.element_width(element_width)) R(clk,A_rows[(no_of_row_by_vector_modules-j)*NI*element_width-1-:element_width*NI],
vector_rows[(no_of_row_by_vector_modules-j)*NI*element_width-1-:element_width*NI],result[(no_of_row_by_vector_modules-j)*element_width-1-:element_width],give_us_all[(no_of_row_by_vector_modules-j-1)],no_of_multiples[(no_of_row_by_vector_modules-j)*32-1-:32],start_row_by_vector[no_of_row_by_vector_modules-j-1],decoder_read_now[no_of_row_by_vector_modules-j-1],reset,you_can_read[no_of_row_by_vector_modules-j-1],I_am_ready[no_of_row_by_vector_modules-j-1]);
	
end
endgenerate


always @(posedge clk)
	begin
	 second_pipeline <= first_pipeline;
	end	 
	
always@(posedge clk)
	begin
		if(give_us_all!=0 && ! (&give_us_all_dash))
			begin
				give_us_all_dash <= give_us_all_dash |	give_us_all;	
			end	
		else if(&give_us_all_dash && !(&give_us_all))	
			begin 
				give_us_all_dash<=0;
			end	
	end	 
	
	always@(posedge clk)
	begin
		if(decoder_read_now!=0 && ! (&decoder_read_now_dash))
			begin
				decoder_read_now_dash <= decoder_read_now_dash | decoder_read_now;	
			end	
		else if(&decoder_read_now_dash && !(&decoder_read_now))	
			begin 
				decoder_read_now_dash<=0;
			end	
	end	




always@(posedge clk) begin
	if(reset)
	begin	
		i<=0;
		first_pipeline <=1;
		start_row_by_vector <=0;
		matrix_by_vector_finished <=0; 
		first_initialization_counter<=0; 
		memories_pre_preprocess<=0;
	end
	
	else if(!start) 		 
		 begin		 
			 i<=0;
			 first_pipeline <=1;  
			 first_initialization_counter<=0;
			 memories_pre_preprocess<=0;
		 end 
	else if(start && first_pipeline && ~first_initialization_counter)
		begin 
			
			memories_pre_preprocess <=1;  
			first_pipeline<=0; 
			first_initialization_counter<=1;
		end	
	else if(start && ~first_pipeline)
		begin  
			memories_pre_preprocess <=0; 
			first_pipeline<=1;
			start_row_by_vector <= -1;
		end		
	 else if(i<(total_with_additional_A/no_of_row_by_vector_modules )&&start && (&give_us_all_dash)  ) 
		
		begin
		memories_pre_preprocess <=1;	
		start_row_by_vector <= -1;
        i<=i+1;	
		if( i ==(total_with_additional_A/no_of_row_by_vector_modules))
			begin	
				 matrix_by_vector_finished <=1;
			end	
		end 
 
	 else /* if((& give_us_all) !=1)  */ // I commented this for issues conecrning the special case
		 begin	
			memories_pre_preprocess <=0; 
			start_row_by_vector <= 0;	 
		 end 
		 
	
	 end	  
	 
	 


always @(posedge clk)
	begin
		if(reset)
			begin
				finish<=0;
				counter2<=0;
			end
			
		else if(!start)
			begin
			finish<=0;
				 
			counter2<=0;	
				
			end
			
			
		else if(counter2==6)
			begin
			finish<=1;	
			
			end
			
		   else if(start)
				counter2<=counter2+1;
			end
		
	
endmodule