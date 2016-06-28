//--------------------------------------------------------------------------------------------
//
// Generated by X-HDL VHDL Translator - Version 2.0.0 Feb. 1, 2011
// ????? ?????? 13 2013 21:06:53
//
//      Input file      : 
//      Component name  : div_nr_wsticky
//      Author          : 
//      Company         : 
//
//      Description     : 
//
//
//--------------------------------------------------------------------------------------------


module div_nr_wsticky(A, B, Q, sticky);
   parameter             NBITS = 6;
   parameter             PBITS = 8;
   input [NBITS-1:0]     A;
   input [NBITS-1:0]     B;
   output [PBITS-1:0]    Q;
   output                sticky;
   
   
   parameter [NBITS-1:0] ZEROS = 0;
   
   reg [NBITS:0]         YY_in[0:PBITS];
   reg [PBITS-1:0]       QQ_in[0:PBITS];
   reg [NBITS:0]         m_cablesIn[0:PBITS];
   wire [NBITS:0]        m_cablesOut[0:PBITS];
   
   wire [PBITS:0]        a_or_s;
   wire [PBITS-1:0]      QQ;
   wire [NBITS-1:0]      YY;
   wire [NBITS-1:0]      XX;
   
   assign XX = A;
   assign YY = B;
   assign Q = QQ;
   
   assign a_or_s[0] = 1'b0;
   always @(*) m_cablesIn[0] <= {1'b0, XX};
   always @(*) YY_in[0] <= {1'b0, YY};
   
   generate
      begin : xhdl0
         genvar                I;
         for (I = 0; I <= PBITS - 1; I = I + 1)
         begin : divisor
            
            a_s #(NBITS) int_mod(.op_a(m_cablesIn[I]), .op_m(YY_in[I]), .as(a_or_s[I]), .outp(m_cablesOut[I]));
         end
      end
   endgenerate
   
   generate
      begin : xhdl1
         genvar                I;
         for (I = 0; I <= PBITS - 1; I = I + 1)
         begin : conex
            assign a_or_s[I + 1] = m_cablesOut[I][NBITS];
            always @(*) m_cablesIn[I + 1] <= {m_cablesOut[I][NBITS - 1:0], 1'b0};
            always @(*) YY_in[I + 1] <= YY_in[I];
            always @(*) QQ_in[I + 1][I] <= a_or_s[I + 1];
            if (I > 0)
            begin : rest
               always @(*) QQ_in[I + 1][I - 1:0] <= QQ_in[I][I - 1:0];
            end
         end
      end
   endgenerate
   
   generate
      begin : xhdl2
         genvar                I;
         for (I = 0; I <= PBITS - 1; I = I + 1)
         begin : quotient
            assign QQ[I] = (~QQ_in[PBITS][PBITS - 1 - I]);
         end
      end
   endgenerate
   
   assign sticky = (m_cablesOut[PBITS - 1] == ZEROS) ? 1'b0 : 
                   1'b1;
   
endmodule