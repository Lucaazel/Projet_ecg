function [ signal,time_axis, N, Fs, R_loc, Q_loc, S_loc, interval_left, interval_right] = QRSDetection()


    clear;
    clc;
    %% Loading signals and freq: ecg and Fs 
    [file,path] = uigetfile('*.mat', 'rt');
    load(fullfile(path, file));
    N = size(ecg,2); % Data length
    time_axis = (1:N)/Fs;


    %% On l'affiche
    plot(time_axis, ecg);
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

    %%five-point differentiation filter (WARNING: Not a causal filter)

    b = [1 2 0 -2 -1];
    a = [ 8/Fs ];
    Y_dec = filter(b, a, Y);

    
    Y_shift = shift(2, Y_dec);  %we shift Y_dec because the filter is not causal.

    %%squaring step:
    M=16;                       %average length of QRS interval (0.08 - 0.10 second)
    s = abs(Y_shift).^2;

    %%moving windows integration (sMWI):
    sMWI = zeros(1, N+2);       %N+2 car il y a un d?callage a cause du filtre acausal
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

    seuil = max(Y_filtre)*0.28 %we chose this threshold, arbitrarily
    figure;
    plot([time_axis 0 0], Y_filtre); 

    %% we build all the intervals (where are the complexes Q, R and S)
    delay = 27;                 %delay of all filters combined
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

    signal = ecg;
    R_loc = RR_indices;
    Q_loc = Q_indices;
    S_loc = S_indices;
    interval_left = left_interval;
    interval_right = right_interval;
    
    

