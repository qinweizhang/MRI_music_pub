% (c) Q. Zhang @ amc 2018

clear; close all; clc
if ispc
    cd('L:\basic\divi\Ima\parrec\jschoormans\Kerry\mp3')
else
    cd('')
end
filename = 'Bach Cello Suite No.1 - Prelude (Yo-Yo Ma).mp3';
[y, Fs ] = audioread(filename);
y = y./max(y(:));
%%
scale = 20;
sample_interval = 1/Fs * 1000; %in ms
GR_wave = y * scale;           %in mT/m
max_str = max(GR_wave,[],1);     %in mT/m
max_sr = max(diff(GR_wave,1)/sample_interval,[], 1);   %in mT/m/ms

disp('==============================================================')
str = sprintf('Sample intervel \t\t: %f ms (%f KHz)',sample_interval, Fs/1000);
disp(str);
str = sprintf('Total duration \t\t\t: %f ms (%f min)',sample_interval * length(y), sample_interval * length(y)/1000/60);
disp(str);
str = sprintf('Scale factor \t\t\t: %f',scale);
disp(str);
str = sprintf('Max GR strength \t\t: %f (x) and %f (y) mT/m',max_str);
disp(str);
str = sprintf('Max GR slewRate \t\t: %f (x) and %f (y) mT/m/ms',max_sr);
disp(str);
disp('==============================================================')



%% filter
plot_flag = 1;
if plot_flag
    figure(1);
    subplot(311);
    plot([1:length(y)]*sample_interval, y); xlabel('ms'); title('original waveform')
end
%fft
spectrum = fft1d(y);

if plot_flag
    figure(2);
    bw = 1/sample_interval; %in kHz
    f_axis = linspace(-bw/2, bw/2, length(spectrum));
    plot(f_axis, abs(spectrum)); xlabel('kHz');title('original spectrum')
end

%hanning window filter
filter_bw = 2; % in kHz (fixed)
str = sprintf('Appling hanning filter: bw: %.2f kHz ...',filter_bw);
disp(str);
passband_idx = find(abs(f_axis)<(filter_bw/2));
% hann_window = hann(length(passband_idx));
hann_window = ones(length(passband_idx),1);
hann_window_matchsize = f_axis * 0;
hann_window_matchsize(passband_idx) = hann_window;
spectrum_filt = bsxfun(@times, spectrum, hann_window_matchsize');
%ifft
y_filterd_long = ifft1d(spectrum_filt);
y_filterd = ifft1d(spectrum_filt(passband_idx,:));
if plot_flag
    hold on
    plot(f_axis, abs(hann_window_matchsize),'g--'); xlabel('kHz');

    figure(3);
    plot(f_axis, abs(spectrum_filt)); xlabel('kHz');
    
    figure(1);
    subplot(312);
    plot([1:length(y)]*sample_interval, y_filterd_long); xlabel('ms'); title('filtered same points #')
    subplot(313);
    plot([1:length(y_filterd)]*(1/filter_bw * 1000), y_filterd); xlabel('ms'); title('filtered: less sample points')
    
end


%%
y_filterd = real(y_filterd);
y_filterd = y_filterd./max(y_filterd(:));
scale = 40; %<<<<<<<<< max gradient strength
sample_interval_filt = 1/(filter_bw * 1000) * 1000; %in ms
slope_dur = 0.5; %in ms - 4 kHz
GR_wave_filt = y_filterd * scale;           %in mT/m
max_str_filt = max(GR_wave_filt,[],1);     %in mT/m
max_sr_filt = max(diff(GR_wave_filt,1)/slope_dur,[], 1);   %in mT/m/ms

disp('==========================After Filtering=================================')
str = sprintf('Sample points \t\t\t: %d ',length(y_filterd));
disp(str);
str = sprintf('Sample intervel \t\t: %f ms (%f KHz)',sample_interval_filt, (filter_bw * 1000)/1000);
disp(str);
str = sprintf('Total duration \t\t\t: %f ms (%f min)',sample_interval_filt * length(y_filterd), sample_interval_filt * length(y_filterd)/1000/60);
disp(str);
str = sprintf('Scale factor \t\t\t: %f',scale);
disp(str);
str = sprintf('Max GR strength \t\t: %f (x) and %f (y) mT/m',max_str_filt);
disp(str);
str = sprintf('Max GR slewRate \t\t: %f (x) and %f (y) mT/m/ms',max_sr_filt);
disp(str);
disp('==============================================================')

%% play filterred: Choise one
% sound(y,Fs);  % stop: clear sound
% sound(real(y_filterd_long(1:end,:)),Fs);  % stop: clear sound
% sound(real(y_filterd), filter_bw * 1000);  % stop: clear sound
% sound(real(GR_wave_filt_export / scale), filter_bw * 1000);  % stop: clear sound
%% export data
% scanner need apx 1s to load 10k lines before startup; don't export too
% much
export_range = 1: min(length(GR_wave_filt), 1.2e6);  %<<<<<set export length here: max 2e6 samples = 5 min
GR_wave_filt_export = GR_wave_filt(export_range,:); 

if length(GR_wave_filt_export)>4e3
    fed_out_weight = repmat(cat(1, ones(length(GR_wave_filt_export)-4e3,1), linspace(1, 0, 4e3)'),[1 2]);   %1s fedout
    GR_wave_filt_export = GR_wave_filt_export .* fed_out_weight;
end


GR_wave_filt_export = cat(1, GR_wave_filt_export, zeros(2e4, 2));  %<<<always append 5s
duration = length(GR_wave_filt_export)*0.25/1000/60; % in min
str = sprintf('%d points to export (duration: %.2f min)', length(GR_wave_filt_export), duration );
disp(str)

GR_wave_filt_3ch = cat(2, GR_wave_filt_export, zeros(length(GR_wave_filt_export),1));
GR_wave_filt_3ch = real(round(GR_wave_filt_3ch * 1e6));  %GR_str_step is 1e-6

[~, fn, ~] = fileparts(filename);
dlmwrite(['L:\basic\divi\Ima\parrec\jschoormans\Kerry\dat\',fn,'[',num2str(duration),' min]_','music_GRwaveforms','_',num2str(filter_bw),'kHz.dat'],GR_wave_filt_3ch,'delimiter','\t');


