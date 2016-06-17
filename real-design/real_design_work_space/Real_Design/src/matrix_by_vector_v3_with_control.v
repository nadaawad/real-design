//-----------------------------------------------------------------------------
//
// Title       : matrix_by_vector2
// Design      : cluster jacoubi
// Author      : Windows User
// Company     : toz fek 
//
//-----------------------------------------------------------------------------
//
// File        : matrix_by_vector2.v
// Generated   : Sun Sep 20 14:38:41 2015
// From        : interface description file
// By          : Itf2Vhdl ver. 1.22
//
//-----------------------------------------------------------------------------
//
// Description : 
//
//-----------------------------------------------------------------------------
`timescale 1 ns / 1 ps

//{{ Section below this comment is automatically maintained
//   and may be overwritten
//{module {matrix_by_vector2}}								 
	
`define zero_filling 32'd0
module matrix_by_vector_v3_with_control (clk,reset,start,mat,vector,out,finish,outsider_read_now);  
	
	
parameter no_of_eqn_per_cluster = 10;
parameter no_of_elements_of_mat = 3*(no_of_eqn_per_cluster-1)+1;
parameter element_width =32;   
parameter no_of_units =4; 
parameter reminder=(no_of_eqn_per_cluster%no_of_units);	
parameter number_of_clusters=1;


parameter NI = 8;
parameter additional = NI-(no_of_eqn_per_cluster%NI); 
parameter total = no_of_eqn_per_cluster+additional ;
         

//}} End of automatically maintained section

// -- Enter your statements here -- //		
input clk,reset,start;	 


wire clk,reset,start;
input mat;
wire[element_width*no_of_elements_of_mat-1:0] mat;  
input vector;
wire[(no_of_eqn_per_cluster)*element_width-1:0]vector;	


wire[element_width*(no_of_units)-1:0] result;
reg sel;

output out;
output reg finish;
wire[2*element_width*no_of_units-1:0] out;

//output reg [(no_of_eqn_per_cluster)*element_width-1:0] out_full;


// DECODER 

decoder_with_control #(.no_of_units(no_of_units)) d (clk,result ,out , (&decoder_read_now),outsider_read_now);

//

integer i=0;
integer counter=0;
integer counter2=0;




reg[no_of_units*(3*element_width)-1:0] input1;
reg[(no_of_units+2)*element_width-1:0] input2;

wire [element_width*((no_of_elements_of_mat+2)+(3*additional))-1:0]modified_mat;
wire[(no_of_eqn_per_cluster+2+additional)*element_width-1:0] modified_vector;



assign modified_mat[element_width*((no_of_elements_of_mat+2)+(3*additional))-1-:element_width*(no_of_elements_of_mat+2)]={`zero_filling,mat,`zero_filling};
assign modified_mat[3*element_width*(total)-1-element_width*(no_of_elements_of_mat+2):0]=0;


assign modified_vector[(no_of_eqn_per_cluster+2+additional)*element_width-1-:(no_of_eqn_per_cluster+2)*element_width]={`zero_filling,vector,`zero_filling};
assign modified_vector[((no_of_eqn_per_cluster+2+additional)*element_width-1)-(no_of_eqn_per_cluster+2)*element_width:0]=0;

//reg [no_of_units-1:0] read_row;
wire [no_of_units-1:0] give_us_all;	 
reg initialization_counter ;  
reg [no_of_units-1:0] number_of_multiples; // if three elements per row , then this is number of multiples of three and equals 1
wire[no_of_units-1:0] decoder_read_now;	 
reg [no_of_units-1:0] start_row_by_vector;	
reg matrix_by_vector_finished ;	 
output wire outsider_read_now;		

genvar j;
generate
for(j=0;j<no_of_units;j=j+1) begin:instantiate_ROW_BY_VECTOR
	
row_by_vector_with_control #(.total(total),.no_of_units(no_of_units)) R(clk,input1[3*element_width*(j+1)-1 -:3*element_width],
input2[element_width*(j+3)-1 -:3*element_width],result[(element_width)*(j+1)-1 -:element_width],give_us_all[j],number_of_multiples[j],start_row_by_vector[j],decoder_read_now[j],reset);
	
end
endgenerate



always@(posedge clk) begin
	if(reset)
	begin	
		i<=0;	
		initialization_counter <=0;
		start_row_by_vector <=0;
		matrix_by_vector_finished <=0;
	end
	
	else if(!start) 		 
		 begin		 
			 i<=0;		 
		 end
	
	 else if(i<(total/no_of_units)&&start && (&give_us_all  || ~initialization_counter)) 
		
		begin
		initialization_counter <=1;
		start_row_by_vector <= -1;	
		
		input1<=modified_mat[element_width*3*(total-no_of_units*i)-1 -:3*element_width*no_of_units] ;
		input2<=modified_vector[element_width*((total+2)-i*(no_of_units))-1 -: (no_of_units+2)*element_width ];
		
		number_of_multiples[0] =1;	
		number_of_multiples[1] =1;
		number_of_multiples[2] =1;
		number_of_multiples[3] =1;	
		
		
        i<=i+1;	
		if( i ==(total/no_of_units))
			begin	
				 matrix_by_vector_finished <=1;
			end	
		end 
 
	 else /* if((& give_us_all) !=1)  */ // I commented this for issues conecrning the special case
		 begin
			start_row_by_vector <= 0;	 
		 end 
		 
	
	 end	  
	 
	 
   
 always@(posedge clk) 
	 begin

		 if(reset)
		 begin
			counter<=0;
			sel<=0;
		 end
		 
	 else if(!start)
		 
		   counter<=0;
   
   else if(!reset&&start)
	
	   begin
		
		   
		 if(counter<7)
			   
			   sel<=1;
			
		   else
			   sel<=~sel;	
	
			   
			   counter<=counter+1;
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