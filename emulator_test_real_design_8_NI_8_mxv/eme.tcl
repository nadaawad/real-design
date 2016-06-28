onerror exit

### START OF SCRIPT ###
set Emulator $::env(Emulator)
puts "Emulator is $Emulator"
project open -prj_file_path veloce.prj
connect $Emulator

hwtrace setmode SSR

configure

targetpower ignore

enableio

hwtrace on

memory download -file b.txt -instance tb.uut.matP_v2.mem
memory download -file b.txt -instance tb.uut.matR.mem
memory download -file b.txt -instance tb.uut.instantiate_P_Emap_8[0].P_Emap_mem.mem
memory download -file b.txt -instance tb.uut.instantiate_P_Emap_8[1].P_Emap_mem.mem
memory download -file b.txt -instance tb.uut.instantiate_P_Emap_8[2].P_Emap_mem.mem
memory download -file b.txt -instance tb.uut.instantiate_P_Emap_8[3].P_Emap_mem.mem

memory download -file a.txt -instance tb.uut.matA.mem
memory download -file col_nos.txt -instance tb.uut.col_nos_memory.mem
memory download -file multiples_matrix.txt -instance tb.uut.multiples_mat.mem

memory download -file memx.v -instance tb.uut.matX.mem

memory download -file parameters.txt -instance tb.uut.parameters_memory.mem







#hwclock setup -clk clk -initval 0 -phase 0 -hi 50 -lo 50

trigger download eme.trigger
trigger position 90


#hwclock run 20  -clk clk;


hwclock run 100  -clk clk;


#hwclock run 100  -clk clk;

clock clicks -milliseconds;
puts "Profiling: Before starting clk"
run
waitfor trigger
clock clicks -milliseconds;
puts "Profiling: After Trigger matured"
#stop

upload -tracedir ./my_first_run.stw -alldata 
memory upload -file out.txt -instance tb.uut.matX.mem

disconnect
exit
