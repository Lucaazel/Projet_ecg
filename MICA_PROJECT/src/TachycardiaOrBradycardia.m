function [ ans ] = TachycardiaOrBradycardia( R_indices, Fs)
%Compute if the patient as a tachycardia or a bradycardia
%   [BPM ans] return the BPM and ans 1: tachycaria ans 0 : Bradycardia ans 2: normal
    N = length(R_indices);
    delta = 0;
    for i=1:N-1
        delta = delta + (1/N)*(R_indices(i+1) - R_indices(i))/Fs;
    end
    delta
    if (delta >= 100)
        ans = [60/delta 1];
    elseif (delta <= 60)
        ans = [60/delta 0];
    else
        ans = [60/delta 2];
    end
end

