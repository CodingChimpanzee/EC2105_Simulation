img_start =imread('C:\Users\강성현\Documents\MATLAB\source_image\start.jpg');
patch = img_start(140:163, 274:297);
patch_t = ones(24,24);
ft_img(24,24,patch,patch_t);
function sum = fourier_trans(u,v,h,w,img)
    sum = 0;
    sum = complex(sum);
    for y = 1:h
        for x= 1:w
          angle = (-2)*pi*(y*u/h+x*v/w);
          ei = complex(img(y,x)*cos(angle),img(y,x)*sin(angle));
          a = real(sum) + real(ei);
          
          disp(142/576);
          %disp(a);
          sum = complex((real(sum) + real(ei))/(h*w),(imag(sum)+imag(ei))/(h*w));
          %disp(ei);
          %disp(sum);
        end
    end
end
function ft_img(h,w,s_img,d_img)
    for u = 1:h
        for v = 1:w
            a= fourier_trans(u-h/2,v-w/2,h,w,s_img);
            d_img(u,v) = a;
        end
    end
end
