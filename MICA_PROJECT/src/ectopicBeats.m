function [ indices, deltas ] = ectopicBeats( R_loc, Fs )
%Determined if the patient has ectopic beats or not
%   return the indices of the ectopic beat(s)
indices = [];
deltas = [];                %All the interval between two consecutive R complex (in s)
for i=1:length(R_loc)-1
    delta = (R_loc(i+1) - R_loc(i))/Fs;
    if (delta > 0.15)
        deltas = [deltas delta];
    else
        deltas = [deltas mean(deltas)];
    end
end

seuil =mean(deltas) + 0.2

for j=1:length(deltas)-1
    if deltas(j) >= seuil
        indices = [indices j];
    end
end


end

