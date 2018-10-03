
function g=Gaussian_filter2d(Filter_size, sigma)
%make a 2d gaussian filter

%size=5; %filter size, odd number
sizex=Filter_size(1);
sizey=Filter_size(1);
g=zeros(sizex,sizey); %2D filter matrix
%gaussian filter
for i=-(sizex-1)/2:(sizex-1)/2
    for j=-(sizey-1)/2:(sizey-1)/2
        x0=(sizex+1)/2; %center
        y0=(sizey+1)/2; %center
        x=i+x0; %row
        y=j+y0; %col
        g(y,x)=exp(-(((x-x0).^2)/(2*sigma(1).^2)...
                     +((y-y0).^2)/(2*sigma(2).^2)...
                     ));   
    end
end
%normalize gaussian filter
sum1=sum(g);
sum2=sum(sum1);
g=g/sum2;
end

%plot 3D
%imagesc(Gaussian_filter2d([100 100],[10 10]))
%set(gca,'dataAspectRatio',[1 1 1])