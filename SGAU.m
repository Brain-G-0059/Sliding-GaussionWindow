function result=SGAU(im, radius, iteration,sigma)
% 基于滑动窗口滤波模板的高斯滤波
%papers: 1) Sub-window Box Filter, Y.Gong, B.Liu, X.Hou, G.Qiu, VCIP2018, Dec.09, Taiwan
%        2) Side Window Filtering, H.Yin, Y.Gong, G.Qiu. CVPR2019
%implemented by Yuanhao Gong

r = radius; %the radius of the side window
k = ones(2*r+1,1)/(2*r+1); %separable kernel
k_L=k; k_L(r+2:end)=0; k_L = k_L/sum(k_L); %half kernel
k_R=flipud(k_L);
m = size(im,1)+2*r; n = size(im,2)+2*r; total = m*n;
[row, col]=ndgrid(1:m,1:n);
offset = row + m*(col-1) - total;
im = double(im);
result = im;
d = zeros(m,n,8,'single');

for ch=1:size(im,3)
    U = padarray(im(:,:,ch),[r,r],'replicate');
    for i = 1:iteration
        %all projection distances
        gausFilter = fspecial('gaussian',[3 3],sigma);
        xy1=[1,1,0;1,1,0;0,0,0];
        xy1=gausFilter.*xy1;
        Y1=xy1./sum(sum(xy1));
        % Y1=mapminmax(xy1,0,1);
        %
        xy2=[0,1,1;0,1,1;0,0,0];
        xy2=gausFilter.*xy2;
        % Y2=mapminmax(xy2,0,1);
        Y2=xy2./sum(sum(xy2));
        
        %
        xy3=[0,0,0;0,1,1;0,1,1];
        xy3=gausFilter.*xy3;
        % Y3=mapminmax(xy3,0,1);
        Y3=xy3./sum(sum(xy3));
        
        %
        xy4=[0,0,0;1,1,0;1,1,0];
        xy4=gausFilter.*xy4;
        % Y4=mapminmax(xy4,0,1);
        Y4=xy4./sum(sum(xy4));
        
        %
        xy5=[1,1,1;1,1,1;0,0,0];
        xy5=gausFilter.*xy5;
        % Y5=mapminmax(xy5,0,1);
        Y5=xy5./sum(sum(xy5));
        
        %
        xy6=[0,0,0;1,1,1;1,1,1];
        xy6=gausFilter.*xy6;
        % Y6=mapminmax(xy6,0,1);
        Y6=xy6./sum(sum(xy6));
        
        %
        xy7=[1,1,0;1,1,0;1,1,0];
        xy7=gausFilter.*xy7;
        % Y7=mapminmax(xy7,0,1);
        Y7=xy7./sum(sum(xy7));
        
        %
        xy8=[0,1,1;0,1,1;0,1,1];
        xy8=gausFilter.*xy8;
        % Y8=mapminmax(xy8,0,1);
        Y8=xy8./sum(sum(xy8));
        d(:,:,8)=imfilter(U,Y8,'replicate')-U;
        d(:,:,7)=imfilter(U,Y7,'replicate')-U;
        d(:,:,6)=imfilter(U,Y6,'replicate')-U;
        d(:,:,5)=imfilter(U,Y5,'replicate')-U;
        d(:,:,4)=imfilter(U,Y4,'replicate')-U;
        d(:,:,3)=imfilter(U,Y3,'replicate')-U;
        d(:,:,2)=imfilter(U,Y2,'replicate')-U;
        d(:,:,1)=imfilter(U,Y1,'replicate')-U;

        
        %         d=abs(d);
        %         figure;
        %         subplot(2,4,1);imshow(d(:,:,1),[]);subplot(2,4,2);imshow(d(:,:,2),[]);
        %         subplot(2,4,3);imshow(d(:,:,3),[]);subplot(2,4,4);imshow(d(:,:,4),[]);
        %         subplot(2,4,5);imshow(d(:,:,5),[]);subplot(2,4,6);imshow(d(:,:,6),[]);
        %         subplot(2,4,7);imshow(d(:,:,7),[]);subplot(2,4,8);imshow(d(:,:,8),[]);
        
        %find the minimal signed distance
        tmp = abs(d);
        [~,ind] = min(tmp,[],3);
        index = offset+total*ind;
        dm = d(index); %signed minimal distance
        
        %update
        U = U + dm;
    end
    result(:,:,ch) = U(r+1:end-r,r+1:end-r);
end
