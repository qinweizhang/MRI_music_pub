function res = fft1d(x)
fctr = size(x,1);
res = zeros(size(x));

% size_x = size(x);
% x = reshape(x, size_x(1), size_x(2), []);
% 
% for ii=1:size(x,3)
%         res(:,:,ii) = sqrt(fctr)*fftshift(ifft2(ifftshift(x(:,:,ii))));
% end
res = (1/sqrt(fctr))*fftshift(fft(ifftshift(x,1),[],1),1);

% res = reshape(res, size_x);

end