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


disconnect
exit
