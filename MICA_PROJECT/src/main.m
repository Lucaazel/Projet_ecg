clear;
clc;
close all;
%% Loading signals and freq: ecg and Fs  
load('/data/ecg_normal_2.mat');
size_ecg = size(ecg);
%plot(ecg);
N = size_ecg(2);

%% On l'affiche
t = linspace(0, 200, N);
%figure(1);

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
M=16;%average length of QRS interval (0.08 - 0.10 second)
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

%% we build all the intervals (where are the complexes Q, R and S)
delay = 27; %delay of all filters combined
interval = [];
i0 = 0;
ifin = 0;
k=0;
RR_indices = [];
while (k<=N)
    k = k+1;
    if (Y_filtre(k) > seuil && i0 == 0) 
       RR_indices = [RR_indices, k-delay]; %We've noticed a delay of 27
       i0 = k-M-delay;
       interval = [interval, i0];
    end
    if (Y_filtre(k) > seuil && Y_filtre(k+1) < seuil)
       ifin = k+M-delay;
       i = ifin;
       interval = [interval, ifin];
       i0=0;
    end
end

%% Q and S detection
size_interval = size(interval);
n_interval = size_interval(2);
Q_indices = [];
S_indices = [];
left_interval = [];
right_interval = [];
%Split the interval into 2 (left and right bounds)
for i=1:n_interval
    if (mod(i, 2) == 1)
        left_interval = [left_interval interval(i)];
    else
        right_interval = [right_interval interval(i)];
    end
end

size_RRindices = size(RR_indices);
n_RRindices = size_RRindices(2);

for i=1:n_RRindices
    q = min(ecg(left_interval(i):RR_indices(i)));
    s = min(ecg(RR_indices(i): right_interval(i)));
    
    for k=left_interval(i):right_interval(i)
        if (ecg(k) == q)
            Q_indices = [Q_indices k];
        elseif (ecg(k) == s)
            S_indices = [S_indices k];
        end
    end
end

%readjustement 

for i=1:length(RR_indices)
    if (ecg(RR_indices(i)+1) > ecg(RR_indices(i)))
        RR_indices(i) = RR_indices(i) + 1;
    elseif ecg(RR_indices(i)-1) > ecg(RR_indices(i))
        RR_indices(i) = RR_indices(i) - 1;
    end
end



%% P and T detection 

% P_indices = [];
% T_indices = [];
% 
% for i=1:n_RRindices 
%     if max((ecg(left_interval(i)-35:RR_indices(i)-10))) < abs(min((ecg(left_interval(i)-35:RR_indices(i)-10))))
%         p = min((ecg(left_interval(i)-35:RR_indices(i)-10)));
%     else 
%         p = max((ecg(left_interval(i)-35:RR_indices(i)-10)));
%     end
%     if max((ecg(RR_indices(i)+10: right_interval(i)+35))) < abs(min((ecg(RR_indices(i)+10: right_interval(i)+35))))
%         t = min((ecg(RR_indices(i)+10: right_interval(i)+35)));
%     else 
%         t = max((ecg(left_interval(i)-35:RR_indices(i)-10)));
%     end
%     
%     for k=left_interval(i)-35:right_interval(i)+35
%         if ((ecg(k)) == p)
%             P_indices = [P_indices k];
%         elseif ((ecg(k)) == t)
%             T_indices = [T_indices k];
%         end
%     end
%     
% end


%% Display


%% Load a signal
[file,path] = uigetfile('*.mat', 'rt');
signal = load(fullfile(path, file));
data = signal.ecg; % Your ecg data
Fs = signal.Fs; % Sampling frequency
N = size(data,2); % Data length
time_axis = (1:N)/Fs;


% figure;
% 
% hold on;
h = plot(time_axis, data); grid on;
TachycardiaOrBradycardia(RR_indices, Fs)
% %plot(time_axis(P_indices),data(P_indices), '*','Color','red'); text(time_axis(P_indices),data(P_indices),' P ','Color','red','FontSize',14);
% plot(time_axis(Q_indices),data(Q_indices), '*','Color','red'); text(time_axis(Q_indices),data(Q_indices),' Q ','Color','red','FontSize',14);
% plot(time_axis(RR_indices),data(RR_indices), '*','Color','red'); text(time_axis(RR_indices),data(RR_indices),' R ','Color','red','FontSize',14);
% plot(time_axis(S_indices),data(S_indices), '*','Color','red'); text(time_axis(S_indices),data(S_indices),' S ','Color','red','FontSize',14);
% %plot(time_axis(T_indices),data(T_indices), '*','Color','red'); text(time_axis(T_indices),data(T_indices),' T ','Color','red','FontSize',14);
% hold off;
% xlabel('Time (s)');
% ylabel('Magnitude');
% title('ECG segment characteristic')
% 




% % % %Display points 
% % % no = 10;
% % % window_ecg = ecg(left_interval(no)-35:right_interval(no)+35);
% % % t_axis = linspace(left_interval(no)-35, right_interval(no)+35, length(window_ecg));
% % % %t_axis = left_interval(1):1/length(window_ecg):right_interval(1);
% % % figure;
% % % hold on;
% % % plot(ecg);
% % % for i=1:n_RRindices
% % %     plot(RR_indices(i), ecg(RR_indices(i)), 'x', 'color', 'red'); text(RR_indices(i), ecg(RR_indices(i)), 'R', 'color', 'red', 'Fontsize', 14);
% % %     plot(Q_indices(i), ecg(Q_indices(i)), 'x', 'color', 'red'); text(Q_indices(i), ecg(Q_indices(i)), 'Q', 'color', 'red', 'Fontsize', 14);
% % %     plot(S_indices(i), ecg(S_indices(i)), 'x', 'color', 'red'); text(S_indices(i), ecg(S_indices(i)), 'S', 'color', 'red', 'Fontsize', 14);
% % % end
% % % 


