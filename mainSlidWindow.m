
clear;
close  all
clc
% I=((imread('D:\Desktop\RTP项目代码20191107\2019111仿真结果\图像3\轮廓检测结果.jpg')));
I=(imread('1.jpg'));
% I=(imread('0.bmp'));
% I=(imread('22.bmp'));


% I=I(1:2:end,1:2:end);

 I = imnoise(I,'salt & pepper',0.1);       % 边缘存在问题
%  I = imnoise(I,'gaussian',0,0.05);        % 边缘没了


im=I;
radius=1;
iteration=10;
% 原始算法
tic
result=SGAU(im, radius, iteration,1);
toc

% %  快速算法
% dif = sum(double(I(:))-result(:));
% str =sprintf('%s%f','原始算法误差为 ',dif);
% disp(str);
s=1;
tic
result0=SideWindowBoxFilter(im, radius, iteration);
toc


% dif = sum(double(I(:))-result0(:));
% str =sprintf('%s%f','快速算法误差为 ',dif);
% disp(str);
% 
figure,imshow([I,result,result0],[]);

% figure,imshow([I,result],[]);

% figure,imshow(double(I)-double(result),[]);

