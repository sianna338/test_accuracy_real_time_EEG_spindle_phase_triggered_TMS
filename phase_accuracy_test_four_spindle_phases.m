%% this tests the accurary of phase dependent stimulation for single and paired pulses

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
Peak = {'S149', 'S145'}
Trough = {'S150', 'S146'}
Rising = {'S151', 'S147'}
Falling = {'S152', 'S148'}


% define trials based on conditions
cfg = [];
cfg.datafile  = 'C:\Users\siann\Downloads\spindle-ppTMS_sub-04_ses-adapt.eeg'
cfg.headerfile = 'C:\Users\siann\Downloads\spindle-ppTMS_sub-04_ses-adapt.vhdr'
cfg.continous = 'yes';
cfg.trialdef.eventtype = 'Stimulus';
cfg.trialdef.eventvalue = horzcat(Peak, Trough, Rising, Falling)
cfg.trialdef.prestim = 0.3; 
cfg.trialdef.poststim = 0.3; 

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

% get phase for each condition 
conditions = {Peak, Trough, Rising, Falling}
% prefixes = {'', '', '', '', '', '', '', ''}; % Adjust prefixes accordingly if needed
trials = cell(1, numel(conditions));
trials_unwrap = cell(1, numel(conditions));
for i = 1:numel(conditions)
    condition_1 = conditions{i}(1)
    condition_2 = conditions{i}(2)
    val_condition_1 = condition_1{1}
    val_condition_2 = condition_2{1}
    trials{i} = xphase(spindle_protocol_test.trialinfo == str2double(val_condition_1(2:end)) ...
        | spindle_protocol_test.trialinfo == str2double(val_condition_2(2:end)), :);
    trials_unwrap{i} = xphase_unwrap(spindle_protocol_test.trialinfo == ...
        str2double(val_condition_1(2:end)) | spindle_protocol_test.trialinfo == ...
        str2double(val_condition_2(2:end)), :);
end

trials_peak = trials{1} 
trials_trough = trials{2}
trials_rising = trials{3} 
trials_falling = trials{4} 

trials_peak_unwrap = trials_unwrap{1} 
trials_trough_unwrap = trials_unwrap{2}
trials_rising_unwrap = trials_unwrap{3} 
trials_falling_unwrap = trials_unwrap{4}

%% plot 
% for not unwrapped phase 
plot_positions = [0.5, 1, 1.5, 2]
colors = {'.g', '.b', '.r', '.k'}
figure;
for i=1:numel(conditions)
    for itrial = 1:size(trials{i},1)
        h(i) = plot(plot_positions(i), trials{i}(itrial, 1495), colors{i}, 'MarkerSize', 10)
        hold on; 
    end 
end 

xlim([0.4 2.1])
ylim([-pi pi]); 
xlabel('targeted phase', 'FontSize', 12); 
ylabel('estimated phase (radians)', 'FontSize', 12); 
title('Estimated Phase vs. Targeted Phase (RMS = 12.6)', 'FontSize', 14); 
h = legend(h, {'Peak', 'Trough', 'Rising', 'Falling'}, ...
           'FontSize', 10, 'Location', 'best'); 
grid on; 


% for unwrapped phase 
plot_positions = [0.5, 1, 1.5, 2]
colors = {'.g', '.b', '.r', '.k'}
figure;
for i=1:numel(conditions)
    for itrial = 1:size(trials{i},1)
        h(i) = plot(plot_positions(i), trials_unwrap{i}(itrial, 1495), colors{i}, 'MarkerSize', 10)
        hold on; 
    end 
end  
hold on; 
t = linspace(0, 2*pi, 1000);
x = (28-24)/2 * cos(t) + (28+24)/2; % create cosine wave for reference 
plot(t/pi + 0.5, x, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]); 
xlim([0.4 2.1])
xlabel('targeted phase', 'FontSize', 12); 
ylabel('estimated phase (radians)', 'FontSize', 12); 
title('Comparison of Estimated Phase (unwrapped) with Cosine Wave, RMS = 12.6', 'FontSize', 14); 
h = legend(h, {'Peak', 'Trough', 'Rising', 'Falling'}, ...
           'FontSize', 10, 'Location', 'best')
grid on; 

%% circular plots
figure;
titles = {'peak stimulation', 'trough stimulation', 'rising flank stimulation', ...
    'falling flank stimulation'}
for i=1:numel(conditions)
    subplot(2,2,i)
    std_peak = circ_std(trials{i}(:, 1495)+pi/2, [], [], 1)
    mean_peak = circ_mean(trials{i}(:, 1495)+pi/2)
    circ_plot(trials{i}(:, 1495)+pi/2, 'pretty', 'bo', true,'linewidth',2,'color','b')
    hold on; 
    text(1.3, 0.8, sprintf('SD: %.2f\nMean: %.2f', std_peak, mean_peak), 'FontSize', 10);
    hold off; 
    title(titles{i}, 'Position',[0.2 1.27])
    text(3, 0.2, ['RMS = 12.6'])
end 
