clc;clear;close all;fclose all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Application 3: 2D FFT and image compression
%% 2D FFT for image compression
img_rgb=imread('cat.png'); % load image
img_gray=rgb2gray(img_rgb); % convert to grayscale

% visualize
figure;
% colored photo
subplot(1,2,1);
imshow(img_rgb);
title('RGB');
% grayscale
subplot(1,2,2);
imshow(img_gray);
title('Grayscale');


%% apply FFT, remove X% of the frequency data with low power
X1=80;
c1=fft2(img_gray);
cutoff=prctile(abs(c1(:)),X1);
cX1=c1;
cX1(abs(c1)<cutoff)=0;


%% inverse FFT, reconstruct image
img_reconstruct=ifft2(cX1);
figure1=figure;
subplot(1,2,1);
imshow(img_gray);
title('Original');
subplot(1,2,2);
imshow(img_reconstruct,[0 255]);
title(['Compressed, discard ',num2str(X1),'% data']);

print(figure1,'-dpng','-r300', ['cat_compressed','.png']);
