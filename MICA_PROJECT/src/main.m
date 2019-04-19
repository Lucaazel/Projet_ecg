clear;
clc;
close all;
%% On charge le signal cardiaque var: ecg et Fs  
load('../data/ecg_normal_4.mat');
size_ecg = size(ecg);
N = size_ecg(2);

%% On l'affiche
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
%%%plot(Y);

%five-point differentiation filter, on a v(n-2) in l'enleve pas mais on sen
%souvient pour la fin 

b = [1 2 0 -2 -1];
a = [ 8/Fs ];
Y_dec = filter(b, a, Y);
Y_shift = decalage(2, Y_dec);

%%%plot(abs(Y_shift));

%%squared:
M=16;
s = abs(Y_shift).^2;

%%sMWI:
sMWI = zeros(1, N+2);
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

%% Avec la m?thode du filtre au lieux de la double somme
h = ones(1, M);
h = 1/M*h;

Y_filtre = filter(h, 1, s);

seuil = max(Y_filtre)*0.3;


signal_qrs = zeros(1, N);

int = [];
R_indice = 0;
for i=1:N
    if (Y_filtre(i) >= seuil)
        int = [i-M, i+M];
        break
    end
end



%On construit les intervalles 
intervalle = [];
i0 = 0;
ifin = 0;
k=0;
while (k<=N)
    k = k+1;
    if (Y_filtre(k) > seuil && i0 == 0)
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

%On stock les indices des R
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
    %bug: on prend un elemenet sur 2 de R
    size_R = size(R);
    R_indice = [];
    for k=1:size_R(2)
        if (mod(k, 2) == 1)
            R_indice = [R_indice, R(k)];
        end
    end
    %
end





