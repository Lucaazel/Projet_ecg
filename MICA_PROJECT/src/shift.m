function [ans] = shift( delay, data )

for i=1:delay
    data = [ 0 data ];
end
ans = data

