all: compile run

HDL_TOP = tb

compile:
	vellib dutwork
	velmap work dutwork
	velanalyze -hdl verilog -vef vefilelist
	velcomp -top work.$(HDL_TOP)

SIM_ARGS = -gui # -gui  # -c (Batch Mode) 

run:
	velrun $(SIM_ARGS) -do eme.tcl

clean:
	rm -rf *work
	rm -rf veloce.wave veloce.med veloce.map  velrunopts.ini tbxbindings.h veloce.log DEBUG_AllBoardsInfos_Host_0.txt CMS.log


