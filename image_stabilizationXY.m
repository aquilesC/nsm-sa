function data=image_stabilizationXY(data)

I=data.I(:,2:end-1,:); %cut the edges

X=[-(size(I,2)-1)/2:(size(I,2)-1)/2];
Y=[-(size(I,1)-1)/2:(size(I,1)-1)/2];
[XX,YY]=meshgrid(X,Y);

Im=(shiftdim(mean(I,1)))';
Xres=(trapz(X,repmat(X,size(Im,1),1).*Im.^2,2)./trapz(X,Im.^2,2));
Im=(shiftdim(mean(I,2),2));
Yres=(trapz(Y,repmat(Y,size(Im,1),1).*Im.^2,2)./trapz(Y,Im.^2,2));
I1=zeros(size(XX,1),size(XX,2),size(I,3));
for it=1:size(I,3)
    I1(:,:,it)=interp2(XX,YY,I(:,:,it),XX+Xres(it)-Xres(1),YY+Yres(it)-Yres(1),'spline');
end
ff=find(isnan(sum(sum(I1,3),1))==0);
I=I1(:,ff,:);

Im=shiftdim(mean(I,2),2);

data.stabil.Xcenter_vibration=mean(abs(diff(Xres)));
 
%%
data.I=I;
data.Im=Im;
data.Xum=data.Xum(1:size(data.I,2));