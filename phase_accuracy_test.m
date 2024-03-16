%% this tests the accurary of phase dependent stimulation for single and paired pulses, yey 

clc 
clear all; 
close all; 

% fieldtrip
path_ft   = 'C:\Users\siann\Downloads\fieldtrip-20231113\fieldtrip-20231113';
addpath(path_ft);
ft_defaults;

sf = 5000; 

%% preparations
% define all the conditions and markers
Peak_single_pulse = 'S149'
Trough_single_pulse = 'S150'
Rising_single_pulse = 'S151'
Falling_single_pulse = 'S152'
Peak_pp_TMS = 'S145'
Trough_pp_TMS = 'S146'
Rising_pp_TMS = 'S147'
Falling_pp_TMS = 'S148'

% define trials based on conditions
cfg = [];
cfg.datafile  = 'C:\Users\siann\Downloads\spindle-ppTMS_sub-04_ses-adapt.eeg'
cfg.headerfile = 'C:\Users\siann\Downloads\spindle-ppTMS_sub-04_ses-adapt.vhdr'
cfg.continous = 'yes';
cfg.trialdef.eventtype      = 'Stimulus';
cfg.trialdef.eventvalue     = {Peak_single_pulse, Trough_single_pulse, Rising_single_pulse, Falling_single_pulse, ...
                                Peak_pp_TMS, Trough_pp_TMS, Rising_pp_TMS, Falling_pp_TMS}; 
cfg.trialdef.prestim        = 1; % idk how much data is needed to extract the phase, 3s too little? 
cfg.trialdef.poststim       = 2; 

cfg = ft_definetrial(cfg);

% read-in data 
cfg.channel = {'C4', 'TP9'}
cfg.reref = 'yes'
cfg.refchannel = 'TP9'
data_eeg = ft_preprocessing(cfg)


% find parts of recording where protocol was tested (8*conditions x 5
% trials, so 40 pulses)
vmrk_file = fileread('C:\Users\siann\Downloads\spindle-ppTMS_sub-04_ses-adapt.vmrk');
vmrk_as_cells = regexp(vmrk_file, '\n', 'split'); % split by line 
mask_start = ~cellfun(@isempty, strfind(vmrk_as_cells, 'test spindle pp TMS'));
index_start = vmrk_as_cells(find(mask_start==1));
line_start = split(index_start(1,1),',',1); 

mask_end = ~cellfun(@isempty, strfind(vmrk_as_cells, 'end test')); 
index_end = vmrk_as_cells(find(mask_end==1));
line_end = split(index_end(1,1),',',1); 

% select the data 
cfg = [];
% cfg.latency = [str2double(line_start(3,1))/sf, str2double(line_end(3,1))/sf]; 
[minValue_start,closestIndex_start] = min(abs(data_eeg.sampleinfo(:,1)-str2double(line_start(3,1))))
[minValue_end,closestIndex_end] = min(abs(data_eeg.sampleinfo(:,1)-str2double(line_end(3,1))))
cfg.trials = [closestIndex_start:closestIndex_end]
cfg.channel = 'C4'
spindle_protocol_test = ft_selectdata(cfg, data_eeg)

% filter the data in the spindle frequency range
% 12-15 Hz, filter order: 4, butterworth
for itrial=1:length(spindle_protocol_test.trial)
    spindle_protocol_test.trial_filtered{itrial} = bandpass(spindle_protocol_test.trial{itrial}(1,:),[12 15],sf)
end 
%% extract the phase for each trial
xh = zeros(length(spindle_protocol_test.trial_filtered),length(spindle_protocol_test.trial_filtered{1}(1,:)))
xphase = zeros(length(spindle_protocol_test.trial_filtered),length(spindle_protocol_test.trial_filtered{1}(1,:)))
xphase_unwrap = zeros(length(spindle_protocol_test.trial_filtered),length(spindle_protocol_test.trial_filtered{1}(1,:)))
for itrial = 1:length(spindle_protocol_test.trial_filtered)
    xh(itrial, :) = hilbert(spindle_protocol_test.trial_filtered{itrial}(1,:));
    xphase_unwrap(itrial, :) = (unwrap(angle(xh(itrial,:))));
    xphase(itrial, :) = angle(xh(itrial,:))
end 

%wrap to -1.5*pi .. 0.5*pi 
% estphase = xphase
% estphase(estphase > 0.5*pi) = estphase(estphase > 0.5*pi)-(2*pi);


trials_peak = xphase(spindle_protocol_test.trialinfo == str2double(Peak_pp_TMS(2:end)) | spindle_protocol_test.trialinfo == str2double(Peak_single_pulse(2:end)), :)
trials_falling = xphase(spindle_protocol_test.trialinfo == str2double(Falling_pp_TMS(2:end)) | spindle_protocol_test.trialinfo == str2double(Falling_single_pulse(2:end)), :)
trials_rising = xphase(spindle_protocol_test.trialinfo == str2double(Rising_pp_TMS(2:end)) | spindle_protocol_test.trialinfo == str2double(Rising_single_pulse(2:end)), :)
trials_trough = xphase(spindle_protocol_test.trialinfo == str2double(Trough_pp_TMS(2:end)) | spindle_protocol_test.trialinfo == str2double(Trough_single_pulse(2:end)), :)

trials_peak_unwrap = xphase_unwrap(spindle_protocol_test.trialinfo == str2double(Peak_pp_TMS(2:end)) | spindle_protocol_test.trialinfo == str2double(Peak_single_pulse(2:end)), :)
trials_falling_unwrap = xphase_unwrap(spindle_protocol_test.trialinfo == str2double(Falling_pp_TMS(2:end)) | spindle_protocol_test.trialinfo == str2double(Falling_single_pulse(2:end)), :)
trials_rising_unwrap = xphase_unwrap(spindle_protocol_test.trialinfo == str2double(Rising_pp_TMS(2:end)) | spindle_protocol_test.trialinfo == str2double(Rising_single_pulse(2:end)), :)
trials_trough_unwrap = xphase_unwrap(spindle_protocol_test.trialinfo == str2double(Trough_pp_TMS(2:end)) | spindle_protocol_test.trialinfo == str2double(Trough_single_pulse(2:end)), :)

% for not unwrapped phase 
figure;
for itrial = 1:size(trials_peak,1)
    h(1) = plot(0.5, trials_peak(itrial, 4995), '.r', 'MarkerSize', 10)
    hold on; 
end 
hold on; 
for itrial = 1:size(trials_falling,1)
    h(2) = plot(1, trials_falling(itrial,4995), '.g', 'MarkerSize', 10)
end 
hold on; 
for itrial = 1:size(trials_trough,1)
    h(3) = plot(1.5, trials_trough(itrial,4995), '.b', 'MarkerSize', 10)
end 
hold on; 
for itrial = 1:size(trials_rising,1)
    h(4) = plot(2, trials_rising(itrial,4995), '.k', 'MarkerSize', 10);
end  
% hold on; 
% t = linspace(0, 2*pi, 1000);
% t = linspace(0, 2*pi, 1000);hist
% x = cos(10*t);
% x = 10*cos(t - pi/2) + 80;
% x = (90-76)/2 * cos(t) + (90+76)/2;
% plot(t/pi + 0.5, x, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]); 
% plot(t*1/2*pi + 0.5, x, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5])
% plot(t, x+80, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
xlim([0.4 2.1])
ylim([-pi pi]); 
xlabel('Time (s)', 'FontSize', 12); 
ylabel('estimated phase (radians)', 'FontSize', 12); 
title('Estimated Phase vs. Targeted Phase (RMS = 12.6', 'FontSize', 14); 
h = legend(h, {'Peak', 'Falling', 'Trough', 'Rising'}, ...
           'FontSize', 10, 'Location', 'best'); 
grid on; 

% in degrees
figure;
for itrial = 1:size(trials_peak,1)
    h(1) = plot(0.5, rad2deg(trials_peak(itrial, 4995))+180, '.r', 'MarkerSize', 10)
    hold on; 
end 
hold on; 
for itrial = 1:size(trials_falling,1)
    h(2) = plot(1, rad2deg(trials_falling(itrial,4995))+180, '.g', 'MarkerSize', 10)
end 
hold on; 
for itrial = 1:size(trials_trough,1)
    h(3) = plot(1.5, rad2deg(trials_trough(itrial,4995))+180, '.b', 'MarkerSize', 10)
end 
hold on; 
for itrial = 1:size(trials_rising,1)
    h(4) = plot(2, rad2deg(trials_rising(itrial,4995))+180, '.k', 'MarkerSize', 10)
end  
% plot(t*1/2*pi + 0.5, x, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5])
% plot(t, x+80, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
xlim([0.4 2.1])
ylim([0 360]); 
xlabel('Time (s)', 'FontSize', 12); 
ylabel('estimated phase (degrees)', 'FontSize', 12); 
title('Estimated Phase vs. Targeted Phase', 'FontSize', 14); 
h = legend(h, {'Peak', 'Falling', 'Trough', 'Rising'}, ...
           'FontSize', 10, 'Location', 'best')
grid on; 

% for unwrapped phase 
figure;
for itrial = 1:size(trials_peak,1)
    h(1) = plot(0.5, trials_peak_unwrap(itrial, 4995), '.r', 'MarkerSize', 10)
    hold on; 
end 
hold on; 
for itrial = 1:size(trials_falling,1)
    h(2) = plot(1, trials_falling_unwrap(itrial,4995), '.g', 'MarkerSize', 10)
end 
hold on; 
for itrial = 1:size(trials_trough,1)
    h(3) = plot(1.5, trials_trough_unwrap(itrial,4995), '.b', 'MarkerSize', 10)
end 
hold on; 
for itrial = 1:size(trials_rising,1)
    h(4) = plot(2, trials_rising_unwrap(itrial,4995), '.k', 'MarkerSize', 10)
end  
hold on; 
% t = linspace(0, 2*pi, 1000);
t = linspace(0, 2*pi, 1000);
% x = cos(10*t);
% x = 10*cos(t - pi/2) + 80;
x = (90-77)/2 * cos(t) + (90+76)/2;
plot(t/pi + 0.5, x, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]); 
% plot(t*1/2*pi + 0.5, x, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5])
% plot(t, x+80, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
xlim([0.4 2.1])
ylim([70 100]); 
xlabel('Time (s)', 'FontSize', 12); 
ylabel('estimated phase (radians)', 'FontSize', 12); 
title('Comparison of Estimated Phase (unwrapped) with Cosine Wave, RMS = 12.6', 'FontSize', 14); 
h = legend(h, {'Peak', 'Falling', 'Trough', 'Rising'}, ...
           'FontSize', 10, 'Location', 'best')
grid on; 
%% circular plots and statistics
figure;

subplot(2,2,1)
std_peak = circ_std(trials_peak(:, 4995)+pi/2, [], [], 1)
mean_peak = circ_mean(trials_peak(:, 4995)+pi/2)
circ_plot(trials_peak(:, 4999)+pi/2, 'pretty', 'bo', true,'linewidth',2,'color','b')
hold on; 
text(1.3, 0.8, sprintf('SD: %.2f\nMean: %.2f', std_peak, mean_peak), 'FontSize', 10);
hold off; 
title('peak stimulation', 'Position',[0.2 1.27])
text(3, 0.2, ['RMS = 12.6'])

subplot(2,2,2)
mean_falling = circ_mean(trials_falling(:, 4995)+pi/2)
std_falling = circ_std(trials_falling(:, 4995), [], [], 1)
circ_plot(trials_falling(:, 4995)+pi/2, 'pretty', 'bo', true,'linewidth',2,'color','b')
hold on;
text(1.6, 0.8, sprintf('SD: %.2f\nMean: %.2f', std_falling, mean_falling), 'FontSize', 10);
hold off;
title('falling stimulation', 'Position',[0.25 1.27])

subplot(2,2,4)
mean_rising = circ_mean(trials_rising(:, 4995)+pi/2)
std_rising = circ_std(trials_rising(:, 4995), [], [], 1)
circ_plot(trials_rising(:, 4995)+pi/2, 'pretty', 'bo', true,'linewidth',2,'color','b')
hold on;
text(1.6, 0.5, sprintf('SD: %.2f\nMean: %.2f', std_rising, mean_rising), 'Color', 'k', 'FontSize', 10);
hold off;
title('rising stimulation', 'Position',[0.251 01.25])

subplot(2,2,3)
mean_trough = circ_mean(trials_trough(:, 4995)+pi/2)
std_trough = circ_std(trials_trough(:, 4995), [], [], 1)
circ_plot(trials_trough(:, 4995)+pi/2, 'pretty', 'bo', true,'linewidth',2,'color','b')
title('trough stimulation', 'Position',[0.25 1.25])
hold on;
text(1.3, 0.5, sprintf('SD: %.2f\nMean: %.2f', std_trough, mean_trough), 'Color', 'k', 'FontSize', 10);
hold off;

%% repeat for RMS = 11.6

cfg = [];
% cfg.latency = [str2double(line_start(3,1))/sf, str2double(line_end(3,1))/sf]; 
line_start = split(index_start(1,2),',',1); 
line_end = split(index_end(1,2),',',1); 
[minValue_start,closestIndex_start] = min(abs(data_eeg.sampleinfo(:,1)-str2double(line_start(3,1))))
[minValue_end,closestIndex_end] = min(abs(data_eeg.sampleinfo(:,1)-str2double(line_end(3,1))))
cfg.trials = [closestIndex_start+1:closestIndex_end]
cfg.channel = 'C4'
spindle_protocol_test_2 = ft_selectdata(cfg, data_eeg)

% filter the data in the spindle frequency range
% 12-15 Hz, filter order: 4, butterworth
spindle_protocol_test_2.trial_filtered = spindle_protocol_test_2.trial
for itrial=1:length(spindle_protocol_test_2.trial)
    spindle_protocol_test_2.trial_filtered{itrial} = bandpass(spindle_protocol_test_2.trial{itrial}(1,:),[12 15],sf)
end 
% spindle_protocol_test_2.trial_filtered = spindle_protocol_test_2.trial_filtered(2:end)
% spindle_protocol_test_2.trialinfo = spindle_protocol_test_2.trialinfo(2:end)
xh = zeros(length(spindle_protocol_test_2.trial_filtered),length(spindle_protocol_test_2.trial_filtered{1}(1,:)))
xphase = zeros(length(spindle_protocol_test_2.trial_filtered),length(spindle_protocol_test_2.trial_filtered{1}(1,:)))
xphase_unwrap = zeros(length(spindle_protocol_test_2.trial_filtered),length(spindle_protocol_test_2.trial_filtered{1}(1,:)))
for itrial = 1:length(spindle_protocol_test_2.trial_filtered)
    xh(itrial, :) = hilbert(spindle_protocol_test_2.trial_filtered{itrial}(1,:));
    xphase_unwrap(itrial, :) = (unwrap(angle(xh(itrial,:))));
    xphase(itrial, :) = angle(xh(itrial,:))
end 
%wrap to -1.5*pi .. 0.5*pi 
% estphase = xphase
% estphase(estphase > 0.5*pi) = estphase(estphase > 0.5*pi)-(2*pi);


trials_peak = xphase(spindle_protocol_test_2.trialinfo == str2double(Peak_pp_TMS(2:end)) | spindle_protocol_test_2.trialinfo == str2double(Peak_single_pulse(2:end)), :)
trials_falling = xphase(spindle_protocol_test_2.trialinfo == str2double(Falling_pp_TMS(2:end)) | spindle_protocol_test_2.trialinfo == str2double(Falling_single_pulse(2:end)), :)
trials_rising = xphase(spindle_protocol_test_2.trialinfo == str2double(Rising_pp_TMS(2:end)) | spindle_protocol_test_2.trialinfo == str2double(Rising_single_pulse(2:end)), :)
trials_trough = xphase(spindle_protocol_test_2.trialinfo == str2double(Trough_pp_TMS(2:end)) | spindle_protocol_test_2.trialinfo == str2double(Trough_single_pulse(2:end)), :)

trials_peak_unwrap = xphase_unwrap(spindle_protocol_test_2.trialinfo == str2double(Peak_pp_TMS(2:end)) | spindle_protocol_test_2.trialinfo == str2double(Peak_single_pulse(2:end)), :)
trials_falling_unwrap = xphase_unwrap(spindle_protocol_test_2.trialinfo == str2double(Falling_pp_TMS(2:end)) | spindle_protocol_test_2.trialinfo == str2double(Falling_single_pulse(2:end)), :)
trials_rising_unwrap = xphase_unwrap(spindle_protocol_test_2.trialinfo == str2double(Rising_pp_TMS(2:end)) | spindle_protocol_test_2.trialinfo == str2double(Rising_single_pulse(2:end)), :)
trials_trough_unwrap = xphase_unwrap(spindle_protocol_test_2.trialinfo == str2double(Trough_pp_TMS(2:end)) | spindle_protocol_test_2.trialinfo == str2double(Trough_single_pulse(2:end)), :)

% for not unwrapped phase 
figure;
for itrial = 1:size(trials_peak,1)
    h(1) = plot(0.5, trials_peak(itrial, 4995), '.r', 'MarkerSize', 10)
    hold on; 
end 
hold on; 
for itrial = 1:size(trials_falling,1)
    h(2) = plot(1, trials_falling(itrial,4995), '.g', 'MarkerSize', 10)
end 
hold on; 
for itrial = 1:size(trials_trough,1)
    h(3) = plot(1.5, trials_trough(itrial,4995), '.b', 'MarkerSize', 10)
end 
hold on; 
for itrial = 1:size(trials_rising,1)
    h(4) = plot(2, trials_rising(itrial,4995), '.k', 'MarkerSize', 10);
end  
% hold on; 
% t = linspace(0, 2*pi, 1000);
% t = linspace(0, 2*pi, 1000);hist
% x = cos(10*t);
% x = 10*cos(t - pi/2) + 80;
% x = (90-76)/2 * cos(t) + (90+76)/2;
% plot(t/pi + 0.5, x, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]); 
% plot(t*1/2*pi + 0.5, x, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5])
% plot(t, x+80, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
xlim([0.4 2.1])
ylim([-pi pi]); 
xlabel('Time (s)', 'FontSize', 12); 
ylabel('estimated phase (radians)', 'FontSize', 12); 
title('Estimated Phase vs. Targeted Phase (RMS = 11.6)', 'FontSize', 14); 
h = legend(h, {'Peak', 'Falling', 'Trough', 'Rising'}, ...
           'FontSize', 10, 'Location', 'best'); 
grid on; 

%% circplot

figure;

subplot(2,2,1)
std_peak = circ_std(trials_peak(:, 4995)+pi/2, [], [], 1)
mean_peak = circ_mean(trials_peak(:, 4995)+pi/2)
circ_plot(trials_peak(:, 4999)+pi/2, 'pretty', 'bo', true,'linewidth',2,'color','b')
hold on; 
text(1.3, 0.8, sprintf('SD: %.2f\nMean: %.2f', std_peak, mean_peak), 'FontSize', 10);
hold off; 
title('peak stimulation', 'Position',[0.1 1.27])
text(3, 0.2, ['RMS = 11.6'])

subplot(2,2,2)
mean_falling = circ_mean(trials_falling(:, 4995)+pi/2)
std_falling = circ_std(trials_falling(:, 4995), [], [], 1)
circ_plot(trials_falling(:, 4995)+pi/2, 'pretty', 'bo', true,'linewidth',2,'color','b')
hold on;
text(1.6, 0.8, sprintf('SD: %.2f\nMean: %.2f', std_falling, mean_falling), 'FontSize', 10);
hold off;
title('falling stimulation', 'Position',[0.1 1.27])

subplot(2,2,4)
mean_rising = circ_mean(trials_rising(:, 4995)+pi/2)
std_rising = circ_std(trials_rising(:, 4995), [], [], 1)
circ_plot(trials_rising(:, 4995)+pi/2, 'pretty', 'bo', true,'linewidth',2,'color','b')
hold on;
text(1.6, 0.5, sprintf('SD: %.2f\nMean: %.2f', std_rising, mean_rising), 'Color', 'k', 'FontSize', 10);
hold off;
title('rising stimulation', 'Position',[0.1 01.25])

subplot(2,2,3)
mean_trough = circ_mean(trials_trough(:, 4995)+pi/2)
std_trough = circ_std(trials_trough(:, 4995), [], [], 1)
circ_plot(trials_trough(:, 4995)+pi/2, 'pretty', 'bo', true,'linewidth',2,'color','b')
title('trough stimulation', 'Position',[0.1 1.25])
hold on;
text(1.3, 0.5, sprintf('SD: %.2f\nMean: %.2f', std_trough, mean_trough), 'Color', 'k', 'FontSize', 10);
hold off;
