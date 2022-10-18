clear; clc, close all;


%% Define array for ECG from anfald3.mat

% load('anfald3.mat');
% ECG = anfald3.data(:,end);
% fs = 200;

%% Define array for ECG from MIT-BIH Arrhythmia database

load('100m.mat');
ECG = val(1,:);
fs = 360;

%%
x_axis = linspace(1,length(ECG),length(ECG)) .* 1/fs;

% Block 1: Bandpass filter

% The lowpass filter:

x = ECG;
y = ECG;
y_lp = ECG;
for n = 1:length(ECG)

    if n - 1 < 1
        y_1 = 0;
    else
        y_1 = y(n - 1);
    end

    if n - 2 < 1
        y_2 = 0;
    else
        y_2 = y(n - 2);
    end
    
    if n - 6 < 1
        x_6 = 0;
    else
        x_6 = x(n - 6);
    end

    if n - 12 < 1
        x_12 = 0;
    else
        x_12 = x(n - 12);
    end

    y_lp(n) = 2 * y_1 - y_2 + 1/32 * (x(n) - 2 * x_6 + x_12);

end


% Verify

% figure
% subplot(2,1,1)
% plot(x_axis, x);
% grid minor
% title("ECG before lowpass filtration")
% xlabel("Time [s]")
% ylabel("Amplitude")
% 
% ax = gca;
% ax.FontSize = 20;
% set(gca,'fontname','times');  % Set it to times
% 
% subplot(2,1,2)
% plot(x_axis, y_lp);
% grid minor
% title("ECG after lowpass filtration")
% xlabel("Time [s]")
% ylabel("Amplitude")
% 
% ax = gca;
% ax.FontSize = 20;
% set(gca,'fontname','times');  % Set it to times

%

% The highpass filter:
x = y_lp;
y_hp = y_lp;
for n = 1:length(ECG)

    if n - 1 < 1
        y_1 = 0;
    else
        y_1 = y_lp(n - 1);
    end

    if n - 16 < 1
        x_16 = 0;
    else
        x_16 = x(n - 16);
    end
    
    if n - 32 < 1
        x_32 = 0;
    else
        x_32 = x(n - 32);
    end

    y_hp(n) = x_16 - 1/32 * ( y_1 + x(n) - x_32 );

end


% Verify
 
% figure
% subplot(2,1,1)
% plot(x_axis, x);
% grid minor
% title("ECG before highpass filtration")
% xlabel("Time [s]")
% ylabel("Amplitude")
% 
% ax = gca;
% ax.FontSize = 20;
% set(gca,'fontname','times');  % Set it to times
% 
% subplot(2,1,2)
% plot(x_axis, y_hp);
% grid minor
% title("ECG after highpass filtration")
% xlabel("Time [s]")
% ylabel("Amplitude")
% 
% ax = gca;
% ax.FontSize = 20;
% set(gca,'fontname','times');  % Set it to times


% Block 2: The derivative operator

x = y_hp;
y_d = y_hp;
for n = 1:length(ECG)

    if n - 1 < 1
        x_1 = 0;
    else 
        x_1 = x(n - 1); 
    end

    if n - 3 < 1
        x_3 = 0;
    else 
        x_3 = x(n - 3); 
    end

    if n - 4 < 1
        x_4 = 0;
    else 
        x_4 = x(n - 4); 
    end

    y_d(n) = 1/8 * ( 2 * x(n) + x_1 - x_3 - 2 * x_4 );

end

% Verify

% figure
% subplot(2,1,1)
% plot(x_axis, x);
% grid minor
% title("ECG before derivative operation")
% xlabel("Time [s]")
% ylabel("Amplitude")
% 
% ax = gca;
% ax.FontSize = 20;
% set(gca,'fontname','times');  % Set it to times
% 
% subplot(2,1,2)
% plot(x_axis, y_d);
% grid minor
% title("ECG after derivative operation")
% xlabel("Time [s]")
% ylabel("Amplitude")
% 
% ax = gca;
% ax.FontSize = 20;
% set(gca,'fontname','times');  % Set it to times


% Block 3: Squaring operation

x = y_d;
y_sq = x.^2;

% Verify

% figure
% subplot(2,1,1)
% plot(x);
% subplot(2,1,2)
% plot(y_sq);
% 
% figure
% subplot(2,1,1)
% plot(x_axis, x);
% grid minor
% title("ECG before squaring")
% xlabel("Time [s]")
% ylabel("Amplitude")
% 
% ax = gca;
% ax.FontSize = 20;
% set(gca,'fontname','times');  % Set it to times
% 
% subplot(2,1,2)
% plot(x_axis, y_sq);
% grid minor
% title("ECG after squaring")
% xlabel("Time [s]")
% ylabel("Amplitude")
% 
% ax = gca;
% ax.FontSize = 20;
% set(gca,'fontname','times');  % Set it to times

% Block 4: Integration with N = 30

x = y_sq;
y_int = y_sq;

w_length = 30 * 1/200;
w_samp = w_length * fs; 

s = w_samp;
for n = 1:length(ECG)

    if n - s < 1
        prior_idx = 1;
    else
        prior_idx = n - s;
    end

    y_int(n) = 1/s * sum( x(prior_idx : n) );

end

%% Verify
% 
figure
subplot(2,1,1)
plot(x_axis, x);
grid minor
title("ECG before integration")
xlabel("Time [s]")
ylabel("Amplitude")

ax = gca;
ax.FontSize = 20;
set(gca,'fontname','times');  % Set it to times

subplot(2,1,2)
plot(x_axis, y_int);
grid minor
title("ECG after integration")
xlabel("Time [s]")
ylabel("Amplitude")

ax = gca;
ax.FontSize = 20;
set(gca,'fontname','times');  % Set it to times


%%

findpeaks(y_int, 'MinPeakHeight', 0, 'MinPeakDistance', 0.2 * fs);

[pks,locs] = findpeaks(y_int, 'MinPeakHeight', 0, 'MinPeakDistance', 0.2 * fs);

%% LEARNING PHASE 

% Learning phase 1 on the first 2s of the ECG peaks to initialize the thresholds

% initializing relevant values
learn_period = 2 * fs;
learn_locs = locs( locs <= learn_period );
learn_pks = pks( 1:length(learn_locs) );

threshold_i1 = mean(learn_pks); 
threshold_i2 =  threshold_i1*0.5;
spki = threshold_i1;
npki = threshold_i2;

qrs = [];
rr_intervals = [];

for k = 1:length(learn_pks)
    % current peak
    peaki = pks(k);

    % 0 for noise (no classification), 1 for normal T1 classification, 2 for searchback T2 classification
    classified_qrs_current_iteration = 0;

    % check if the current location is 200ms after the last detected QRS
    % (refractory period)
    timer = 0;
    if length(qrs) > 0 
        timer = locs(k) - qrs(end);
    else
        timer = round(fs * 0.2) + 1;
    end

    if peaki > threshold_i1 && timer > round( fs * 0.2)
        % potential QRS detected
        classified_qrs_current_iteration = 1;

        % save the peak location
        qrs = [qrs, locs(k)];
        
        % save rr intervals after two succesfully QRS values were detected
        if length(qrs) > 2
            rr_intervals = [rr_intervals, (qrs(end) - qrs(end-1))];
        end

    end
           
    % check if the detected QRS in this iteration is a T-wave
    % if rr interval less then 360 ms
    if ~isempty(rr_intervals) && length(qrs) > 2

        if rr_intervals(end) < round(0.36*fs) && classified_qrs_current_iteration > 0 && length(qrs) > 1

            % finding the maximum slopes in relevant intervals

            % interval: (location of the previously detected QRS + refractory period) -> (current location of the detected QRS)
            slope_current = max(abs(diff(y_int(qrs(end-1)+round(0.2*fs):qrs(end)))));

            % interval: (location of the detected QRS before the previously detected QRS + refractory period) -> (location of the previously detected QRS)
            slope_last = max(abs(diff(y_int(qrs(end-2)+round(0.2*fs):qrs(end-1))))); 
            
            % t wave check
            if slope_current < slope_last/2
                % delete the already logged QRS values of this iteration
                % from the arrays
                qrs(end) = [];
                rr_intervals(end) = [];

                % no real QRS detected -> noise peak
                classified_qrs_current_iteration = 0;
            end
        end
    end

    % update the noise and signal peak depending on if a qrs was detected
    if classified_qrs_current_iteration == 0
        % noise peak
        npki = 0.125 * peaki + 0.875 * npki;
    elseif classified_qrs_current_iteration == 1
        % signal peak
        spki = 0.125 * peaki + 0.875 * spki;
    end

    % update thresholds
    threshold_i1 = npki + 0.25 * (spki - npki);
    threshold_i2 = 0.5 * threshold_i1;
end


%% REAL ITERATION (after learning phase) 
threshold_i1_array = [];
threshold_i2_array = [];

% reset arrays and intialize rr-thresholds
qrs = [];
rr_intervals = [];
rr_average2 = 0;
rr_average1 = 0;
rr_missed_limit = 0;

for k = 1:length(pks)

    % save current threshold to array for the plots
    threshold_i1_array = [threshold_i1_array, threshold_i1];
    threshold_i2_array = [threshold_i2_array, threshold_i2];

    % initializing and updating rr interval 1 and 2
    if ~isempty(rr_intervals)

        % rr_average 1
        if length(rr_intervals) >= 8
            rr_n = rr_intervals(end - 7: end);
            rr_average1 = mean(rr_n);
        else
            rr_n = rr_intervals;
            rr_average1 = mean(rr_n);
        end
   
        % initializing rr_average2 for the limits below
        if rr_average2 == 0
            rr_average2 = rr_intervals(end);
        end

        % updating limits
        rr_low_limit = 0.92 * rr_average2;
        rr_high_limit = 1.16 * rr_average2;

        % updating rr_average after limits have been intialized
        all_rr_in_threshold_array = rr_intervals(find(rr_intervals >= rr_low_limit & rr_intervals <= rr_high_limit));
        if ~isempty(all_rr_in_threshold_array)
            if length(all_rr_in_threshold_array) >= 8
                rr_n_in_threshold = all_rr_in_threshold_array(end - 7: end);
                rr_average2 = mean(rr_n_in_threshold);
            else 
                rr_n_in_threshold = all_rr_in_threshold_array;
                rr_average2 = mean(rr_n_in_threshold);
            end
        end
        
        % checking if heartbeat is regular or irregular
        if all( rr_n >= rr_low_limit & rr_average1 <= rr_high_limit )
            % regular
            rr_average2 = rr_average1;
        else
            % irregular
            threshold_i1 = threshold_i1 * 0.5;
        end
        
        % updating missed limit threshold
        rr_missed_limit = 1.66 * rr_average2;
    end
    
    % current peak
    peaki = pks(k);

    % 0 for noise (no classification), 1 for normal T1 classification, 2 for searchback T2 classification
    classified_qrs_current_iteration = 0;
    
    % check if the current location is 200ms after the last detected QRS
    % (refractory period)
    timer = 0;
    if length(qrs) > 0
        timer = locs(k) - qrs(end);
    else
        timer = round(fs * 0.2) + 1;
    end


    if peaki > threshold_i1 && timer > round(fs * 0.2)

        % potential QRS detected
        classified_qrs_current_iteration = 1;

        % save the potential qrs location to the qrs array 
        qrs = [qrs, locs(k)];
        
        % save the rr interval of the most recent beat to the rr array
        if length(qrs) > 2
            rr_intervals = [rr_intervals, qrs(end) - qrs(end-1)];
        end

    else
        % check if searchback needs to be applied

        % placeholder variable to update the potential signal peak in line
        % XX
        search_back_max_peak = 0;

        % rr missed limit has to be initialized
        if rr_missed_limit ~= 0

            % interval between last detected qrs and current iteration
            last_interval = locs(k) - qrs(end);

            % check if this interval is above the rr_missed_limit threshold
            if rr_missed_limit < last_interval
                % initialize search back

                % time to which we want to search back
                peak_gap_timer = qrs(end) + 0.2 * fs;

                % time from which we want to search back
                current_location = locs(k);

                % array to save all the relevant peak locations in this interval
                relevant_peak_indexes = [];

                % save all the peak locations in this interval to an array
                for ll = 1:length(locs)
                    if (locs(ll) >= peak_gap_timer && locs(ll) <= current_location)
                        relevant_peak_indexes = [relevant_peak_indexes, ll];
                    end
                end
                
                % get the max peak value in this interval and its location
                search_back_max_peak = max(pks(relevant_peak_indexes));
                search_back_max_peak_t = locs(find(pks == search_back_max_peak));

                % check if this max peak value is above the threshold I2
                if search_back_max_peak > threshold_i2 && locs(k) - qrs(end) > 0.2 * fs

                    % potential QRS detected in search back
                    classified_qrs_current_iteration = 2;

                    % save the potential qrs location to the qrs array
                    qrs = [qrs, search_back_max_peak_t];

                    % update rr intervals
                    if length(qrs) > 2
                        rr_intervals = [rr_intervals, (qrs(end) - qrs(end-1))];
                    end
                end
            end
        end
    end

    % check if the detected QRS in this iteration is a T-wave
    % if rr interval less then 360 ms
    if ~isempty(rr_intervals) && length(qrs) > 2

        if rr_intervals(end) < round(0.36*fs) && classified_qrs_current_iteration > 0 && length(qrs) > 1

            % finding the maximum slopes in relevant intervals

            % interval: (location of the previously detected QRS + refractory period) -> (current location of the detected QRS)
            slope_current = max(abs(diff(y_int(qrs(end-1)+round(0.2*fs):qrs(end)))));

            % interval: (location of the detected QRS before the previously detected QRS + refractory period) -> (location of the previously detected QRS)
            slope_last = max(abs(diff(y_int(qrs(end-2)+round(0.2*fs):qrs(end-1))))); 
            
            % t wave check
            if slope_current < slope_last/2

                % delete the already logged QRS values of this iteration
                % from the arrays
                qrs(end) = [];
                rr_intervals(end) = [];

                % no real QRS detected -> noise peak
                classified_qrs_current_iteration = 0;
            end
        end
    end

    % update the noise and signal peak depending on if a qrs was detected 
    if classified_qrs_current_iteration == 0
        % noise peak
        npki = 0.125 * peaki + 0.875 * npki;
    elseif classified_qrs_current_iteration == 1
        % signal peak
        spki = 0.125 * peaki + 0.875 * spki;
    elseif classified_qrs_current_iteration == 2 && search_back_max_peak ~= 0
        % searchback signal peak
        spki = 0.25 * search_back_max_peak + 0.75 * spki;
    end

    % update thresholds
    threshold_i1 = npki + 0.25 * (spki - npki);
    threshold_i2 = 0.5 * threshold_i1;
end



%% Verify the total Pan-Thompkins algorithm
length(qrs)

idx_axis = 0:1:length(y_int) - 1;

figure
hold on
plot(idx_axis * 1/fs,  y_int, 'linewidth', 1);

plot( idx_axis(qrs) * 1/fs, y_int(qrs), 'o', 'linewidth', 2);

plot( locs * 1/fs, threshold_i1_array, '--', 'linewidth', 3, 'color', 'g');
plot( locs * 1/fs, threshold_i2_array, '--', 'linewidth', 3, 'color', '#EDB120');
xlim([0, length(y_int) * 1/fs])

legend('','Found peaks', 'Threshold 1', 'Threshold 2', 'Location','northwest')
xlim([1230, 1310]);
grid minor

%title("Pan-Tompkins algorithm on the given ECG after integration")

title("Pan-Tompkins algorithm on MIT-BIH Arrhythmia (No. 222) after integration")
xlabel("Time [s]")
ylabel("Amplitude")

ax = gca;
ax.FontSize = 20;
set(gca,'fontname','times');  % Set it to times





