function [ ans ] = TachycardiaOrBradycardia( R_indices, Fs)
%Compute if the patient as a tachycardia or a bradycardia
%   [BPM ans] return the BPM and ans 1: tachycaria ans 0 : Bradycardia ans 2: normal
    N = length(R_indices);
    delta = 0;
    for i=1:N-1
        delta = delta + (1/N)*(R_indices(i+1) - R_indices(i))/Fs;
    end
    bpm = round(60/delta);
    if (bpm >= 100)
        ans = [bpm 1];
    elseif (bpm <= 60)
        ans = [bpm 0];
    else
        ans = [bpm 2];
    end
end

