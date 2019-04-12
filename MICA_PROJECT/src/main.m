clear;
clc;
%%On charge le signal cardiaque var: ecg et Fs  
load('/data/ecg_normal_1.mat');
size_ecg = size(ecg);
N = size_ecg(2);

%%On l'affiche
t = linspace(0, 200, N);
%figure(1);
%plot(t, ecg);

%%Filtres


b_low_pass = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a_low_pass = [1 -2 1];

b_high_pass = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a_high_pass = [1 -1];

fvtool(b_low_pass, a_low_pass);
fvtool(b_high_pass, a_high_pass);

X_low_pass = filter(b_low_pass, a_low_pass, ecg);
Y = filter(b_high_pass, a_high_pass, X_low_pass);
plot(Y);

%five-point differentiation filter, on a v(n-2) in l'enleve pas mais on sen
%souvient pour la fin 

b = [1 2 0 -2 -1];
a = [ 8/Fs ];
Y_dec = filter(b, a, Y);
%Y_shift = delayseq(Ydec, 2);

plot(abs(Y_dec));

%%squared:
M=10
s = abs(Y_dec).^2;

%%sMWI:
sMWI = zeros(1, N);
su = 0;
%for n=1:N-M+1
%    su(n) = ;
%end


