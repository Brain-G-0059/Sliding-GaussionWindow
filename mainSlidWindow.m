
clear;
close  all
clc
% I=((imread('D:\Desktop\RTP��Ŀ����20191107\2019111������\ͼ��3\���������.jpg')));
I=(imread('1.jpg'));
% I=(imread('0.bmp'));
% I=(imread('22.bmp'));


% I=I(1:2:end,1:2:end);

 I = imnoise(I,'salt & pepper',0.1);       % ��Ե��������
%  I = imnoise(I,'gaussian',0,0.05);        % ��Եû��


im=I;
radius=1;
iteration=10;
% ԭʼ�㷨
tic
result=SGAU(im, radius, iteration,1);
toc

% %  �����㷨
% dif = sum(double(I(:))-result(:));
% str =sprintf('%s%f','ԭʼ�㷨���Ϊ ',dif);
% disp(str);
s=1;
tic
result0=SideWindowBoxFilter(im, radius, iteration);
toc


% dif = sum(double(I(:))-result0(:));
% str =sprintf('%s%f','�����㷨���Ϊ ',dif);
% disp(str);
% 
figure,imshow([I,result,result0],[]);

% figure,imshow([I,result],[]);

% figure,imshow(double(I)-double(result),[]);

