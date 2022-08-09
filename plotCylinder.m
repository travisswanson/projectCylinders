function plotCylinder(x,y,z,dx,dy,dz,radius,c2p)
noCirPoints = 180;

tCir = linspace(0,2*pi(),noCirPoints)';

XOrtho = cos(tCir);
YOrtho = sin(tCir);
%unit vector at base of cylinder, pointing up cylinder axis
aV = [dx,dy,dz]./(sqrt(dx.^2+dy.^2+dz.^2));


%function to find orthgonal vectors
oVz = @(v,x,y)((-v(1).*x - v(2).*y)./v(3));



%calculate set of orthgonal vectors
ZOrtho = oVz(aV,XOrtho,YOrtho);

%unit-ify the orthgonal vectors
uov = [XOrtho,YOrtho,ZOrtho]./sqrt(XOrtho.^2 + YOrtho.^2 + ZOrtho.^2);

  
%donot re unit-fy, you only want the horizontal component, not the
%renormalized horizontal component

%using only the X and Y components, find circle coods in plane of
%interest
xaC = x(1) + uov(:,1).*radius;
yaC = y(1) + uov(:,2).*radius;
zaC = z(1) + uov(:,3).*radius;

xbC = x(2) + uov(:,1).*radius;
ybC = y(2) + uov(:,2).*radius;
zbC = z(2) + uov(:,3).*radius;

surf([xaC,xbC],[yaC,ybC],[zaC,zbC],'FaceColor',c2p,'EdgeColor',c2p)

%

end