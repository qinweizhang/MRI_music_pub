%%
clear sound
p = y(3e5:end,:);
p(:,2) = 0;
x_idx = 1:length(p);

%% decompose channels
% ttt(:,1) = sin(x_idx'*pi/1000).*p(:,1); 
% ttt(:,2) = cos(x_idx'*pi/1000).*p(:,1); 
% %% split channels
% mask = mod(round(x_idx/1e6),2);
% ttt(:,1) = p(:,1) .* mask';
% ttt(:,2) = p(:,1) .* (1 - mask);
%% simple offset/scale
ttt = p+1;
%% single pitch
% ttt = 10 * ones(length(p),2);

%%
sound(ttt,Fs); 

figure(1);
subplot(221); plot(p); title('true waveform') 
subplot(222); plot(abs(fft1d(p))); title('spectrum')
subplot(223); plot(ttt); title('current waveform') 
subplot(224); plot(abs(fft1d(ttt))); title('current spectrum')