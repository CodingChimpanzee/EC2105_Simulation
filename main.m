%-------------------------------------------------------------------------%
% Fourier Transform

startimg = imread('C:\Users\강성현\Documents\MATLAB\source_image\start.jpg');
startgray = rgb2gray(startimg);
xtart = 270;
ytart = 132;
w = 39;
h = 38;
face = imcrop(startgray, [xtart ytart w-1 h-1]);    % Make image Gray

result = fourier_transform(face, w, h);
result = result/(w*h);
% Get the result patch
faceout = (log(abs(result)));
subplot(1,2,1), imagesc(faceout);
% Print the result of face fourier transform

total = fourier_transform(startgray, 480, 360);
total = total/(360*480);
totalout = (log(abs(total)));
subplot(1,2,2), imagesc(totalout);
% Print the result of start.jpg fourier transform

%-------------------------------------------------------------------------%
% Inverse Fourier Transform

face_recov = inverse_fourier(result, w, h);
face_recov = (abs(face_recov));
subplot(2,2,3), imshow(face_recov);
% Print the recovered face

total_recov = inverse_fourier(total, 480, 360);
total_recov = (abs(total_recov));
subplot(2,2,4), imshow(total_recov);
% Print the recovered start.jpg

%-------------------------------------------------------------------------%
% Phase Correlation

value = zeros(480, 360);
%for a = 1 : (480-w)
%    for b = 1 : (360-h)
%        patch = imcrop(startgray, [a b w-1 h-1]);
%        G = fourier_transform(patch, w, h);
%        R = (result.*conj(G))./abs(result.*conj(G));
%        value(a, b) = sum(abs(inverse_fourier(R, w, h)), 'all');
%    end
%end




