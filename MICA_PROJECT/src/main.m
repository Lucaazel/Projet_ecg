clear;
clc;
close all;
%% Loading signals and freq: ecg and Fs  
load('../data/ecg_normal_4.mat');
size_ecg = size(ecg);
N = size_ecg(2);

%% On l'affiche
t = linspace(0, 200, N);
%figure(1);
plot(ecg);

%%Filters


b_low_pass = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a_low_pass = [1 -2 1];

b_high_pass = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a_high_pass = [1 -1];

%fvtool(b_low_pass, a_low_pass);
%fvtool(b_high_pass, a_high_pass);

X_low_pass = filter(b_low_pass, a_low_pass, ecg);
Y = filter(b_high_pass, a_high_pass, X_low_pass);
%%%plot(Y);

%five-point differentiation filter (WARNING: Not a causal filter)

b = [1 2 0 -2 -1];
a = [ 8/Fs ];
%fvtool(b, a);
Y_dec = filter(b, a, Y);

%we shift Y_dec because the filter is not causal.
Y_shift = shift(2, Y_dec);

%%%plot(abs(Y_shift));

%%squaring step:
M=16;
s = abs(Y_shift).^2;

%%moving windows integration (sMWI):
sMWI = zeros(1, N+2);%N+2 car il y a un d?callage a cause du filtre acausal
for n=1:N+2
    subtotal = 0;
    for i=0:M-1
        if (n-i) > 0
            subtotal = subtotal + s(n-i);
        else
            subtotal = 0;
        end
    end 
    sMWI(n) = (1/M)*subtotal;
end

%% It is easier like this with the filter() function: 
h = ones(1, M);
h = 1/M*h;

Y_filtre = filter(h, 1, s);

seuil = max(Y_filtre)*0.32; %we chose this threshold, arbitrarily

%% we build all the intervals (where are the complexes P, Q, R and S)
intervalle = [];
i0 = 0;
ifin = 0;
k=0;
RR_indices = [];
while (k<=N)
    k = k+1;
    if (Y_filtre(k) > seuil && i0 == 0) 
       RR_indices = [RR_indices, k-27]; %We've noticed a delay of 27
       i0 = k-M;
       intervalle = [intervalle, i0];
    end
    if (Y_filtre(k) > seuil && Y_filtre(k+1) < seuil)
       ifin = k+M;
       i = ifin;
       intervalle = [intervalle, ifin];
       i0=0;
    end
end

%% R indices are saved: -------Useless now thanks to RR_indices
R = [];
size_intervalle = size(intervalle);
indice = 0;
for i=-1:size_intervalle(2)-3
    debut = intervalle(i+2);
    fin = intervalle(i+3);
    ecg_int = ecg(:, debut-23:fin-23);
    max_ecg_int = max(abs(ecg_int));
    
    indice = find(abs(ecg)==max_ecg_int);
    R = [R, indice];
    %bug: we take one element out of 2
    size_R = size(R);
    R_indice = [];
    for k=1:size_R(2)
        if (mod(k, 2) == 1)
            R_indice = [R_indice, R(k)];
        end
    end
    %
end





