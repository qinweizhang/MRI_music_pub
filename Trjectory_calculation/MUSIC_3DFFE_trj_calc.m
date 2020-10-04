% function to calc k-space trjactory of the muisc GR waveform
% 
% (C) Kerry zhang  2018 @AMC

%% Settings
clc; close all;

fn = 'Bach_yoyoma_2kHz.dat';

%==========const.=================
%Units are: ms, mT, m
gamma       = 42.576;   %in /(mT * ms) or (KHz/mT)
gr_dwell    = 0.5;   %in ms  2KHz
gr_str_step = 1e-6;
TR = 8; %ms
dummy_TR = 20; % every dummy_TR+1 TR one AQ. 
AQ_start = 0.9871; %ms
AQ_end = 0.9871+5.0425; %ms
FA = 10; %degree

gr_samples_per_TR =  TR / gr_dwell * (dummy_TR+1);
disp([num2str(gr_samples_per_TR),' smaples per TR']);

max_gr_length = 480000; % <500,000 this matters for CDAS
%% Load 1ch GR Waveform (audio waveform)
% gr is the ready-to-use waveform, paired with GR_str_step of 1e-6
% but only for sound effect. no trajectory design
if(ispc)
    gr = dlmread(['L:\basic\divi\Ima\parrec\jschoormans\Kerry\dat\',fn],'\t');
else
    gr = dlmread(['/home/qzhang/lood_storage/divi/Ima/parrec/jschoormans/Kerry/dat/',fn],'\t');
end
gr_str = gr .* gr_str_step;


%% decompose gr into three directions for K-space sampling: ku
single_ch_gr = gr(:,1);
gr_3ch = zeros(length(single_ch_gr),3);
% Options:
%   2ch_radial,             1ch_const,          1ch_to_m, , 
%   1ch_to_3ch_kushball,    1ch_to_2ch_spiral,  2ch_untouched,
%   1ch_to_3ch_cone,        1ch_to_2ch_radial,  1ch_sin,
%   EPI_X                   EPI_2Dradial
option = 'EPI_2Dradial'; 

switch option
    case '1ch_to_3ch_kushball'
        % change GR direction once per TR
        disp('decomposing GR using the Kushball option')
        %----------initilization-------------------------
        TR_total = 2000; %asume 2000 TR used (5min with TR = 12ms)
        k_trj = zeros(TR_total, gr_samples_per_TR, 3);

        for idx = 1:length(single_ch_gr)
           n_TR = floor(idx/gr_samples_per_TR);
           %set polar angles
           theta = n_TR * 16.95229 / 180 * pi; % golden angle
           phi = n_TR / TR_total * (-pi)  + pi / 2;  % pi/2 --> -pi/2
           %clac compositions
           gr_3ch(idx, 1) = cos(phi) * cos(theta) * single_ch_gr(idx);  %Gx
           gr_3ch(idx, 2) = cos(phi) * sin(theta) * single_ch_gr(idx);  %Gy
           gr_3ch(idx, 3) = sin(phi) * single_ch_gr(idx);  %Gz
           %calc K-space trajecotry
           if(n_TR + 1 <= TR_total)
               local_sample_n = mod(idx, gr_samples_per_TR) + 1;
               if (local_sample_n + 1 <= gr_samples_per_TR)
                   for dir = 1:3
                     k_trj(n_TR+1, local_sample_n + 1, dir) =  k_trj(n_TR+1, local_sample_n, dir)  + ...
                         gr_3ch(idx, dir) .* gr_str_step .* gamma .* gr_dwell;  
                   end
               end
           end
        end
    case '1ch_to_3ch_cone'
        %TODO: continousely change GR direction 
        disp('decomposing GR using the continous 3D spiral option')
    case '1ch_to_2ch_spiral'
        disp('option 1ch_to_2ch_spiral used')
        for idx = 1:length(single_ch_gr)
           gr_3ch(idx, 1) = cos(idx/200) * single_ch_gr(idx);  %Gx
           gr_3ch(idx, 2) = sin(idx/200) * single_ch_gr(idx);  %Gy
        end
    case '1ch_to_2ch_radial'
         % change GR direction once per TR
        disp('decomposing GR using the 2D radial option')
        %----------initilization-------------------------
        TR_total = 2000; %asume 2000 TR used (5min with TR = 12ms)
        k_trj = zeros(TR_total, gr_samples_per_TR, 3);

        for idx = 1:length(single_ch_gr)
           n_TR = floor(idx/gr_samples_per_TR);
           %set polar angles
           theta = n_TR * 16.95229 / 180 * pi; % golden angle
           phi = 0;
           %clac compositions
           gr_3ch(idx, 1) = cos(phi) * cos(theta) * single_ch_gr(idx);  %Gx
           gr_3ch(idx, 2) = cos(phi) * sin(theta) * single_ch_gr(idx);  %Gy
           gr_3ch(idx, 3) = sin(phi) * single_ch_gr(idx);  %Gz
           %calc K-space trajecotry
           if(n_TR + 1 <= TR_total)
               local_sample_n = mod(idx, gr_samples_per_TR) + 1;
               if (local_sample_n + 1 <= gr_samples_per_TR)
                   for dir = 1:3
                     k_trj(n_TR+1, local_sample_n + 1, dir) =  k_trj(n_TR+1, local_sample_n, dir)  + ...
                         gr_3ch(idx, dir) .* gr_str_step .* gamma .* gr_dwell;  
                   end
               end
           end
        end
    case '1ch_to_m'
        disp('option 1ch_to_m used')
        gr_3ch(:,1) = single_ch_gr;
    case '1ch_sin'
        str = 3; %mT/m
        cycle = 19; %number of samples
        mgg = sprintf('option 1ch_sin, cycle: %.2f samples', cycle);
        disp(mgg);
        gr_3ch(:,1) = str / gr_str_step * sin([1:length(gr_3ch)]' / cycle * 2 * pi);  %use constent readout to test recon pipeline
    case '1ch_const'
        str = 3; %mT/m
        max_stew_rate = 20; %mT/m/ms
        max_slope_step = 0.25 * max_stew_rate; %mT/m
        gr_3ch(:,1) = str / gr_str_step * ones(size(gr_3ch(:,1)));  %use constent readout to test recon pipeline
        slope = 0:max_slope_step:str;
        gr_3ch(1:length(slope),1) = slope / gr_str_step;
        
    case '2ch_radial'
        disp('option 2ch_radial used')
        str = 1; %mT/m
        max_stew_rate = 20; %mT/m/ms
        max_slope_step = 0.25 * max_stew_rate; %mT/m
        single_ch_gr = str / gr_str_step * ones(size(gr_3ch(:,1)));  %use constent readout to test recon pipeline
        slope = 0:max_slope_step:str;
        single_ch_gr(1:length(slope)) = slope / gr_str_step;
        
        % change GR direction once per TR
        disp('decomposing GR using the 2D radial option')
        %----------initilization-------------------------
        TR_total = 2000; %asume 2000 TR used (5min with TR = 12ms)
        k_trj = zeros(TR_total, gr_samples_per_TR, 3);

        for idx = 1:length(single_ch_gr)
           n_TR = floor(idx/gr_samples_per_TR);
           %set polar angles
           theta = n_TR * 16.95229 / 180 * pi; % golden angle
           phi = 0;
           %clac compositions
           gr_3ch(idx, 1) = cos(phi) * cos(theta) * single_ch_gr(idx);  %Gx
           gr_3ch(idx, 2) = cos(phi) * sin(theta) * single_ch_gr(idx);  %Gy
           gr_3ch(idx, 3) = sin(phi) * single_ch_gr(idx);  %Gz
           %calc K-space trajecotry
           if(n_TR + 1 <= TR_total)
               local_sample_n = mod(idx, gr_samples_per_TR) + 1;
               if (local_sample_n + 1 <= gr_samples_per_TR)
                   for dir = 1:3
                     k_trj(n_TR+1, local_sample_n + 1, dir) =  k_trj(n_TR+1, local_sample_n, dir)  + ...
                         gr_3ch(idx, dir) .* gr_str_step .* gamma .* gr_dwell;  
                   end
               end
           end
        end
    case 'EPI_X'
        disp('option EPI_X used');
        %         segment = [0 -25 -25 0 25 50 50 50 50 25 0 -25 -25 -25 0 0]'/20 ;
        segment = [0  0  -1.5 0 3 -3  3 -3 3 -3  3 -3  3 -3   3 -3 ]'/gr_str_step * 5;
%         segment = [0 -1.5 3 -3  3 -3 3 -3  3 -3  3 -3   3 -3  3 -3]'/gr_str_step * 5;
        temp = repmat(segment, [ceil(length(gr_3ch)/length(segment)),1]);
        gr_3ch(:, 1) = temp(1:length(gr_3ch)); 
        
        trj_seg = segment * 0;
        for ii = 2:length(trj_seg)
            trj_seg(ii) =  trj_seg(ii-1) + gamma * gr_3ch(ii-1) * gr_str_step * gr_dwell; % in 1/m
            %                           /(mT * ms)      mT/m                       ms   
        end
        AQ_window = [AQ_start AQ_end];  %AQ start-end time point(ms)
        figure(11);
        subplot(211); plot([0:(length(segment)-1)]*gr_dwell, segment*gr_str_step); hold on; plot(AQ_window, 0*AQ_window, 'r-','LineWidth',2); title('GR waveform segment (input)'); xlabel('ms'); ylabel('GRx mT/m')
        subplot(212); plot([0:(length(segment)-1)]*gr_dwell, trj_seg); title('Trajectory segment (input)'); hold on; plot(AQ_window, 0*AQ_window, 'r-','LineWidth',2); xlabel('ms'); ylabel('kx 1/m')
     case 'EPI_2Dradial'
        disp('option EPI_2Dradial used');
        %         segment = [0 -25 -25 0 25 50 50 50 50 25 0 -25 -25 -25 0 0]'/20 ;
%         segment = [0  0 0 0 -1.5 0 3 -3  3 -3 3 -3  3 -3  3 -3 ]'/gr_str_step * 5;
        %         segment = [0 -1.5 3 -3  3 -3 3 -3  3 -3  3 -3   3 -3  3 -3]'/gr_str_step * 5;
        segment = [0  0 0 0  -1.5 0 3 3 0 -3 -3 0 3 3 0 -3]'/gr_str_step * 5;
        single_ch_gr = repmat(segment, [ceil(length(gr_3ch)/length(segment)),1]);
        
        
        %----------initilization-------------------------
        TR_total = 2000; %asume 2000 TR used; only for display purposes
        k_trj = zeros(TR_total, gr_samples_per_TR, 3);
        for idx = 1:length(single_ch_gr)
            n_TR = floor(idx/(gr_samples_per_TR)); %consider dummy TR
            %set polar angles
            theta = n_TR * 16.95229 / 180 * pi; % golden angle
            phi = 0;
            %clac compositions
            gr_3ch(idx, 1) = cos(phi) * cos(theta) * single_ch_gr(idx);  %Gx
            gr_3ch(idx, 2) = cos(phi) * sin(theta) * single_ch_gr(idx);  %Gy
            gr_3ch(idx, 3) = sin(phi) * single_ch_gr(idx);  %Gz
            %calc K-space trajecotry
            if(n_TR + 1 <= TR_total)
                local_sample_n = mod(idx, gr_samples_per_TR) + 1;
                if (local_sample_n + 1 <= gr_samples_per_TR)
                    for dir = 1:3
                        k_trj(n_TR+1, local_sample_n + 1, dir) =  k_trj(n_TR+1, local_sample_n, dir)  + ...
                            gr_3ch(idx, dir) .* gr_str_step .* gamma .* gr_dwell;
                    end
                end
            end
        end
        
        AQ_window = [AQ_start AQ_end];  %AQ start-end time point(ms)
        figure(11);
        plot_tr_n = 50;
        subplot(411); plot([0:(gr_samples_per_TR-1)]*gr_dwell, gr_3ch([1:gr_samples_per_TR]+plot_tr_n*gr_samples_per_TR,1)*gr_str_step); hold on; plot(AQ_window, 0*AQ_window, 'r-','LineWidth',2); title('GR x waveform segment (input)'); xlabel('ms'); ylabel('GRx mT/m')
        subplot(412); plot([0:(gr_samples_per_TR-1)]*gr_dwell, gr_3ch([1:gr_samples_per_TR]+plot_tr_n*gr_samples_per_TR,2)*gr_str_step); hold on; plot(AQ_window, 0*AQ_window, 'r-','LineWidth',2); title('GR y waveform segment (input)'); xlabel('ms'); ylabel('GRx mT/m')
        subplot(413); plot([0:(gr_samples_per_TR-1)]*gr_dwell, squeeze(k_trj(plot_tr_n,:,1))); title('Trajectory x segment (input)'); hold on; plot(AQ_window, 0*AQ_window, 'r-','LineWidth',2); xlabel('ms'); ylabel('kx 1/m')
        subplot(414); plot([0:(gr_samples_per_TR-1)]*gr_dwell, squeeze(k_trj(plot_tr_n,:,2))); title('Trajectory y segment (input)'); hold on; plot(AQ_window, 0*AQ_window, 'r-','LineWidth',2); xlabel('ms'); ylabel('kx 1/m')

    otherwise
        disp('option otherwise used')
        gr_3ch = gr;
end


%% plot
figure(1);
plot(gr_3ch); title('gradient waveform')
switch option
    case '1ch_to_3ch_kushball'
        figure(3); 
        plot3(0,0,0,'ro','MarkerFaceColor','r'); hold on; grid on;
        xlim([-1 1]); ylim([-1 1]);zlim([-1 1]);
        for tr = 1:2000
        plot3(squeeze(k_trj(tr,[1, 60],1)), squeeze(k_trj(tr,[1, 60],2)),squeeze(k_trj(tr,[1, 60],3)));
        % plot3([min(k_trj(tr,:,1)), max(k_trj(tr,:,1))], squeeze(k_trj(tr,[1, 60],2)),squeeze(k_trj(tr,[1, 60],3)));

        end
        % plot_range= 1:length(k_all);
        % scatter(k_all(plot_range,1),k_all(plot_range,2), [],  abs(gr_str(plot_range,1))); 
        % xlabel('kx (1/m)'); ylabel('ky (1/m)')
        % % plot3(k_all(plot_range,1),k_all(plot_range,2),k_all(plot_range,3),'o');
    otherwise
        
end


%% exporting the results
% scanner need apx 1s to load 10k lines before startup; don't export too
% much
if length(gr_3ch)>max_gr_length
    gr_3ch  = gr_3ch(1:max_gr_length,:);
end

duration = length(gr_3ch)*0.25/1000/60; % in min

[~, fnn, ~] = fileparts(fn);

if(ispc)
    dlmwrite(['L:\basic\divi\Ima\parrec\jschoormans\Kerry\dat\FFE_',option,'__',fnn,'[',num2str(duration),' min]_','music_GRwaveforms.dat'],gr_3ch,'delimiter','\t');
else
%     dlmwrite(['/home/qzhang/lood_storage/divi/Users/qzhang/qzhang_user/MRI_music/dat/FFE_',fnn,'[',num2str(duration),' min]_','music_GRwaveforms.dat'],gr_3ch,'delimiter','\t');
    dlmwrite(['/home/qzhang/lood_storage/divi/Ima/parrec/jschoormans/Kerry/dat/FFE_',option,'__',fnn,'[',num2str(duration),' min]_','music_GRwaveforms.dat'],gr_3ch,'delimiter','\t');
    
end


