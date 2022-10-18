clear; clc;
load("anfald3.mat")

fs = 200;
f_cutoff = [5, 15]/(fs/2);
n = 6;
ftype = 'bandpass';
[b,a] = butter(n, f_cutoff, ftype);

fieldnames = anfald3.channels();
%fvtool(b,a)

% design frequeny range
f = linspace(0, 200, 1301);
legend_names = [];

figure(); hold on
for jj = 1:21
    current_seizure_part = anfald3.data(3300:end,jj);
    filtered_current_seizure_part = filtfilt(b, a, current_seizure_part);
    filtered_current_seizure_part_transform = fft(filtered_current_seizure_part);
    filtered_current_seizure_part_transform_abs = abs(filtered_current_seizure_part_transform);

    plot(f, filtered_current_seizure_part_transform_abs);
end
xlabel('Frequency (Hz)')
ylabel('Amplitude (A.U.)')
xlim([0, 25])
legend(fieldnames)
title('Frequency domain for each EEG signal with bandpass filter')
hold off

fh = findall(0,'Type','Figure');
txt_obj = findall(fh,'Type','text');
set(txt_obj,'fontname','times','FontSize', 15);  % Set it to times
