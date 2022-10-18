%clear; clc;
set(gcf,'renderer','Painters')
%% Loading the data
nn1=importdata("nn1.xlsx");
nn2=importdata("nn2.xlsx");
nn3=importdata("nn3.xlsx");
name_list = ["TimeStamp", "LifetouchHeartRate", "LifetouchRespirationRate", "OximeterSpO2","OximeterPulse", "BloodPressureSystolic", "BloodPressureDiastolic", "BloodPressurePulse" ];

%% Create a patient structs
% Create a struct 1
nn1_struct = struct();
for ii = 1:length(nn1.textdata)
    fieldname = name_list(ii);
    rowvalues = nn1.data(:,ii);
    nn1_struct.(fieldname{1}) = rowvalues;
end

% Create a struct 2
nn2_struct = struct();
for ii = 1:length(nn2.textdata)
    fieldname = name_list(ii);
    rowvalues = nn2.data(:,ii);
    nn2_struct.(fieldname{1}) = rowvalues;
end

% Create a struct 3
nn3_struct = struct();
for ii = 1:length(nn3.textdata)
    fieldname = name_list(ii);
    rowvalues = nn3.data(:,ii);
    nn3_struct.(fieldname{1}) = rowvalues;
end

%% cleaning data 
% high values
nn1_struct.LifetouchHeartRate = filter_high(nn1_struct.LifetouchHeartRate);
nn2_struct.LifetouchHeartRate = filter_high(nn2_struct.LifetouchHeartRate);
nn3_struct.LifetouchHeartRate = filter_high(nn3_struct.LifetouchHeartRate);

nn1_struct.LifetouchRespirationRate = filter_high(nn1_struct.LifetouchRespirationRate);
nn2_struct.LifetouchRespirationRate = filter_high(nn2_struct.LifetouchRespirationRate);
nn3_struct.LifetouchRespirationRate = filter_high(nn3_struct.LifetouchRespirationRate);

% null values
nn1_struct.OximeterPulse = filter_null(nn1_struct.OximeterPulse);
nn2_struct.OximeterPulse = filter_null(nn2_struct.OximeterPulse);
nn3_struct.OximeterPulse = filter_null(nn3_struct.OximeterPulse);

nn1_struct.OximeterSpO2 = filter_null(nn1_struct.OximeterSpO2);
nn2_struct.OximeterSpO2 = filter_null(nn2_struct.OximeterSpO2);
nn3_struct.OximeterSpO2 = filter_null(nn3_struct.OximeterSpO2);


%% plotting data
plot_patient(nn1_struct, "1")
plot_patient(nn2_struct, "2")
plot_patient(nn3_struct, "3")

fh = findall(0,'Type','Figure');
txt_obj = findall(fh,'Type','text');
set(txt_obj,'fontname','times','FontSize', 15);  % Set it to times

%% create frequency table
frequency_table_patient_1 = create_frequency_table(nn1_struct);
frequency_table_patient_2 = create_frequency_table(nn2_struct);
frequency_table_patient_3 = create_frequency_table(nn3_struct);

percentage_frequency_table_patient_1 = compute_percentages(frequency_table_patient_1);
percentage_frequency_table_patient_2 = compute_percentages(frequency_table_patient_2);
percentage_frequency_table_patient_3 = compute_percentages(frequency_table_patient_3);

function percent_table = compute_percentages(table)
table_values = table(:,2:end);
for i=1:4
  table_values.(i) = str2double(table_values{:,i});
end

for i=1:5
    rowsum = sum(table_values{i,:});
    for j=1:4
        table_values{i, j} = table_values{i, j}/rowsum;
    end
end
percent_table = table_values;
end


%% reusable functions

% high filter function > 1000
function array_out = filter_high(array_in)
upper_filter = 1000;
array_out = [];
for ii = 1:length(array_in)
    current_value = array_in(ii) ;
    if (current_value < upper_filter)
        array_out(ii) = current_value;
    else
        array_out(ii) = NaN;
    end
end
array_out = array_out.';
end

% null filter function == 0
function array_out = filter_null(array_in)
null_filter = 0;
array_out = [];
for ii = 1:length(array_in)
    current_value = array_in(ii);
    if (current_value >= null_filter)
        array_out(ii) = current_value;
    else
        array_out(ii) = NaN;
    end
end
array_out = array_out.';
end

% function counting the scores for heart rate
function array_out = heart_rate_counter(hr_array_in)
score_array = [0, 0, 0, 0];
for i = 1:length(hr_array_in)
    if hr_array_in(i) < 90.5 && hr_array_in(i) >= 50.5
        score_array(1) = score_array(1)+1;
    end
    if hr_array_in(i) < 50.5 && hr_array_in(i) >= 40.5
        score_array(2) = score_array(2)+1;
    end
     if hr_array_in(i) < 110.5 && hr_array_in(i) >= 90.5
        score_array(2) = score_array(2)+1;
    end
    if hr_array_in(i) < 130.5 && hr_array_in(i) >= 110.5
        score_array(3) = score_array(3)+1;
    end
    if hr_array_in(i) < 40.5 || hr_array_in(i) >= 130.5
        score_array(4) = score_array(4)+1;
    end
array_out = score_array;
end
end

% function counting the scores for respiratory rate
function array_out = respiratory_rate_counter(hr_array_in)
score_array = [0, 0, 0, 0];
for i = 1:length(hr_array_in)
    if hr_array_in(i) < 20.5 && hr_array_in(i) >= 11.5
        score_array(1) = score_array(1)+1;
    end
    if hr_array_in(i) < 11.5 && hr_array_in(i) >= 8.5
        score_array(2) = score_array(2)+1;
    end
    if hr_array_in(i) < 24.5 && hr_array_in(i) >= 20.5
        score_array(3) = score_array(3)+1;
    end
    if hr_array_in(i) < 8.5 || hr_array_in(i) >= 24.5
        score_array(4) = score_array(4)+1;
    end
array_out = score_array;
end
end

% function counting the scores for oxygen saturation
function array_out = oxygen_saturation_counter(hr_array_in)
score_array = [0, 0, 0, 0];
for i = 1:length(hr_array_in)
    if hr_array_in(i) >= 95.5
        score_array(1) = score_array(1)+1;
    end
    if hr_array_in(i) < 95.5 && hr_array_in(i) >= 93.5
        score_array(2) = score_array(2)+1;
    end
    if hr_array_in(i) < 93.5 && hr_array_in(i) >= 91.5
        score_array(3) = score_array(3)+1;
    end
    if hr_array_in(i) < 91.5
        score_array(4) = score_array(4)+1;
    end
array_out = score_array;
end
end

% function counting the scores for systolic blood pressure
function array_out = systolic_blood_pressure_counter(hr_array_in)
score_array = [0, 0, 0, 0];
for i = 1:length(hr_array_in)
    if hr_array_in(i) < 219.5 && hr_array_in(i) >= 110.5
        score_array(1) = score_array(1)+1;
    end
    if hr_array_in(i) < 110.5 && hr_array_in(i) >= 100.5
        score_array(2) = score_array(2)+1;
    end
    if hr_array_in(i) < 100.5 && hr_array_in(i) >= 90.5
        score_array(3) = score_array(3)+1;
    end
    if hr_array_in(i) < 90.5 || hr_array_in(i) >= 219.5
        score_array(4) = score_array(4)+1;
    end
array_out = score_array;
end
end

% create and export score table into LaTeX with n-values and %
function table_out = create_frequency_table(patient_struct)
headers = {'Vital Sign','Score 0','Score 1','Score 2','Score 3'};
empty_table = cell2table(cell(0,5),'VariableNames', headers);
table = [empty_table; num2cell(["LifetouchHeartRate", heart_rate_counter(patient_struct.LifetouchHeartRate);])];
table = [table; num2cell(["OximeterPulse", heart_rate_counter(patient_struct.OximeterPulse);])];
table = [table; num2cell(["LifetouchRespirationRate", respiratory_rate_counter(patient_struct.LifetouchRespirationRate);])];
table = [table; num2cell(["OximeterSpO2", oxygen_saturation_counter(patient_struct.OximeterSpO2);])];
table = [table; num2cell(["BloodPressureSystolic", systolic_blood_pressure_counter(patient_struct.BloodPressureSystolic);])];
table_out = table;
check_table_value_length(patient_struct.LifetouchHeartRate, heart_rate_counter(patient_struct.LifetouchHeartRate));
check_table_value_length(patient_struct.OximeterPulse, heart_rate_counter(patient_struct.OximeterPulse));
check_table_value_length(patient_struct.LifetouchRespirationRate, respiratory_rate_counter(patient_struct.LifetouchRespirationRate));
check_table_value_length(patient_struct.OximeterSpO2, oxygen_saturation_counter(patient_struct.OximeterSpO2));
check_table_value_length(patient_struct.BloodPressureSystolic, systolic_blood_pressure_counter(patient_struct.BloodPressureSystolic));
end

function check_table_value_length(data_array, table_array)
    if length(data_array(~isnan(data_array))) ~= sum(table_array)
        error = "ERROR - Table data count does not reflect the real data"
    end 
end


function plot_patient(patient_struct, patient_id)
%% data plot of Lifetouch Heart Rate
% https://de.mathworks.com/help/matlab/ref/yline.html
% https://de.mathworks.com/matlabcentral/answers/454747-changing-color-of-the-matlab-plot-at-different-ranges
figure;
x = linspace(1, 1, length(patient_struct.LifetouchHeartRate)+1);
X=[0:length(patient_struct.LifetouchHeartRate),fliplr(0:length(patient_struct.LifetouchHeartRate))];
% score 0
y1=50.5*(x);  
y2=90.5*(x);
Y_0=[y1,fliplr(y2)]; 
fill(X,Y_0,'g');
hold on
% score 1
y3=40.5*(x);
y4=50.5*(x);
y5=90.5*(x);
y6=110.5*(x);
Y_1_1=[y3,fliplr(y4)]; 
Y_1_2=[y5,fliplr(y6)];
fill(X,Y_1_1,'y');
hold on
fill(X,Y_1_2,'y');
hold on
% score 2
y7=110.5*(x);  
y8=130.5*(x);
Y_2 = [y7,fliplr(y8)]; 
fill(X,Y_2,[0.8500 0.3250 0.0980]);
hold on
%score 3
y9=30*(x);
y10=40.5*(x);
y11=130.5*(x);
y12=180*(x);
Y_3_1 = [y9,fliplr(y10)]; 
Y_3_2 = [y11,fliplr(y12)]; 
fill(X,Y_3_1,'r');
hold on
fill(X,Y_3_2,'r');
hold on
alpha(0.2)
grid
plot(1:length(patient_struct.LifetouchHeartRate),patient_struct.LifetouchHeartRate, "color", "k")
hold off
xlabel('Time [min]')
ylabel('Heart Rate [bpm]')
title(join(['Lifetouch Heart Rate (patient: ',patient_id, ')']))
legend('score 0 - [51-90]','score 1 - [41-50] or [91-110]', '', 'score 2 - [111-130]', 'score 3 - [<41 or >130]')
xlim([0, length(patient_struct.LifetouchHeartRate)])
ylim([30, 180])
set(gcf,'renderer','Painters')

%% data plot of Oximeter pulse
figure;
x = linspace(1, 1, length(patient_struct.OximeterPulse)+1);
X=[0:length(patient_struct.OximeterPulse),fliplr(0:length(patient_struct.OximeterPulse))];
% score 0
y1=50.5*(x);  
y2=90.5*(x);
Y_0=[y1,fliplr(y2)]; 
fill(X,Y_0,'g');
hold on
% score 1
y3=40.5*(x);
y4=50.5*(x);
y5=90.5*(x);
y6=110.5*(x);
Y_1_1=[y3,fliplr(y4)]; 
Y_1_2=[y5,fliplr(y6)];
fill(X,Y_1_1,'y');
hold on
fill(X,Y_1_2,'y');
hold on
% score 2
y7=110.5*(x);  
y8=130.5*(x);
Y_2 = [y7,fliplr(y8)]; 
fill(X,Y_2,[0.8500 0.3250 0.0980]);
hold on
%score 3
y9=20*(x);
y10=40.5*(x);
y11=130.5*(x);
y12=250*(x);
Y_3_1 = [y9,fliplr(y10)]; 
Y_3_2 = [y11,fliplr(y12)];
fill(X,Y_3_1,'r');
hold on
fill(X,Y_3_2,'r');
hold on
alpha(0.2)
grid
plot(1:length(patient_struct.OximeterPulse),patient_struct.OximeterPulse, 'color', 'k')
xlabel('Time [min]')
ylabel('Pulse [bpm]')
title(join(['Oximeter Pulse (patient: ',patient_id, ')']))
legend('score 0 - [51-90]','score 1 - [41-50] or [91-110]', '', 'score 2 - [111-130]', 'score 3 - [<41 or >130]')
xlim([0, length(patient_struct.OximeterPulse)])
ylim([20, 250])
set(gcf,'renderer','Painters')

%% data plot of Lifetouch Respiratory Rate
figure;
x = linspace(1, 1, length(patient_struct.LifetouchRespirationRate)+1);
X=[0:length(patient_struct.LifetouchRespirationRate),fliplr(0:length(patient_struct.LifetouchRespirationRate))];
% score 0
y1=11.5*(x);  
y2=20.5*(x);
Y_0=[y1,fliplr(y2)]; 
fill(X,Y_0,'g');
hold on
% score 1
y3=8.5*(x);
y4=11.5*(x);
Y_1_1=[y3,fliplr(y4)]; 
fill(X,Y_1_1,'y');
hold on
% score 2
y7=20.5*(x);  
y8=24.5*(x);
Y_2 = [y7,fliplr(y8)]; 
fill(X,Y_2,[0.8500 0.3250 0.0980]);
hold on
%score 3
y9=5*(x);
y10=8.5*(x);
Y_3_1 = [y9,fliplr(y10)];
y11=24.5*(x);
y12=35*(x);
Y_3_2 = [y11,fliplr(y12)]; 
fill(X,Y_3_1,'r');
hold on
fill(X,Y_3_2,'r');
hold on
alpha(0.2)
grid
plot(1:length(patient_struct.LifetouchRespirationRate),patient_struct.LifetouchRespirationRate, "color", "k")
xlabel('Time [min]')
ylabel('Respiration Rate [breaths per minute]')
title(join(['Lifetouch Respiration Rate (patient: ',patient_id, ')']))
legend('score 0 - [12-20]','score 1 - [9-11]', 'score 2 - [21-24]', 'score 3 - [<9 or >24]')
xlim([0, length(patient_struct.LifetouchRespirationRate)])
ylim([5, 35])
set(gcf,'renderer','Painters')

%% data plot of Oximeter Sp02
figure;
x = linspace(1, 1, length(patient_struct.OximeterSpO2)+1);
X=[0:length(patient_struct.OximeterSpO2),fliplr(0:length(patient_struct.OximeterSpO2))];
% score 0
y1=95.5*(x);  
y2=100*(x);
Y_0=[y1,fliplr(y2)]; 
fill(X,Y_0,'g');
hold on
% score 1
y3=93.5*(x);
y4=95.5*(x);
Y_1_1=[y3,fliplr(y4)]; 
fill(X,Y_1_1,'y');
hold on
% score 2
y7=91.5*(x);  
y8=93.5*(x);
Y_2 = [y7,fliplr(y8)]; 
fill(X,Y_2,[0.8500 0.3250 0.0980]);
hold on
%score 3
y9=50*(x);
y10=91.5*(x);
Y_3_1 = [y9,fliplr(y10)]; 
fill(X,Y_3_1,'r');
hold on
alpha(0.2)
grid
plot(1:length(patient_struct.OximeterSpO2),patient_struct.OximeterSpO2, "color", "k")
xlabel('Time [min]')
ylabel('SpO2 [%]')
title(join(['Oximeter SpO2 (patient: ',patient_id, ')']))
legend('score 0 - [>95]','score 1 - [94-95]', 'score 2 - [92-93]', 'score 3 - [<92]')
xlim([0, length(patient_struct.OximeterSpO2)])
ylim([50, 100])
set(gcf,'renderer','Painters')


%% data plot of Bloodpressure Systolic
% finding the indices that is not a NaN (would not plot otherwise (in this
% case))
idx=find(patient_struct.BloodPressureSystolic>0);
figure;
x = linspace(1, 1, length(patient_struct.OximeterSpO2)+1);
X=[0:length(patient_struct.BloodPressureSystolic),fliplr(0:length(patient_struct.BloodPressureSystolic))];
% score 0
y1=110.5*(x);  
y2=219.5*(x);
Y_0=[y1,fliplr(y2)]; 
fill(X,Y_0,'g');
hold on
% score 1
y3=100.5*(x);
y4=110.5*(x);
Y_1_1=[y3,fliplr(y4)]; 
fill(X,Y_1_1,'y');
hold on
% score 2
y7=90.5*(x);  
y8=100.5*(x);
Y_2 = [y7,fliplr(y8)]; 
fill(X,Y_2,[0.8500 0.3250 0.0980]);
hold on
%score 3
y9=70*(x);
y10=90.5*(x);
Y_3_1 = [y9,fliplr(y10)]; 
y11=219.5*(x);
y12=250*(x);
Y_3_2 = [y11,fliplr(y12)]; 
fill(X,Y_3_1,'r');
fill(X,Y_3_2,'r');
hold on
alpha(0.2)
grid
plot(idx,patient_struct.BloodPressureSystolic(idx), "color", "k")
xlabel( 'Time [min]')
ylabel('Systolic Bloodpressure [mmHg]')
title(join(['Bloodpressure Systolic (patient: ',patient_id, ')']))
legend('score 0 - [111-219]','score 1 - [101-110]', 'score 2 - [91-100]', '', 'score 3 - [<91 or >219]')
xlim([0, length(patient_struct.BloodPressureSystolic)])
ylim([70, 250])
set(gcf,'renderer','Painters')
end




