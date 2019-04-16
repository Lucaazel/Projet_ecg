function [rep] = decalage( retard, data )

for i=1:retard
    data = [ 0 data ];
end
rep = data

