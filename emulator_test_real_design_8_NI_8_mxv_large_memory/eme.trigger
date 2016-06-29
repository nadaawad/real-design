domain domain0 @(posedge clk);

state init : initial, store
begin
  next ( tb.halt == 1);
end

state ok: trigger,stop;


