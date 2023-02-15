function T = projectCylinders(T)

    warning('off','MATLAB:polyshape:repairedBySimplify')

    x = T.x;
    y = T.y;
    z = T.z;
    radius = T.radius;
    noCylinders = T.noCylinders;
    cLength = T.cLength;
    dx = T.dx;
    dy = T.dy;
    dz = T.dz;
    theta = T.theta;
    BO = T.BO;
    maxBO =T.maxBO;


    % Projected area of cylinders and shadow fraction
    
    noCirPoints = 360;
    
    tCir = linspace(0,2*pi(),noCirPoints)';
    
    XOrtho = cos(tCir);
    YOrtho = sin(tCir);
    
    % [XOrtho, YOrtho]= poly2cw(XOrtho,YOrtho);
    
    %unit vector at base of cylinder, pointing up cylinder axis
    aV = [dx,dy,dz]./(sqrt(dx.^2+dy.^2+dz.^2));
    bV = -aV; %unit vector looking down from top circle (but not translated)
    
    %function to find orthgonal vectors
    oVz = @(v,x,y)((-v(1).*x - v(2).*y)./v(3));
    
    minZ = zeros(1,noCylinders);
    maxZ = zeros(1,noCylinders);
    
    %for each cylinder
    for idx = 1:noCylinders
    
        %in the case there's no horizontal movement of the cylinder ends, it's
        %area is a circle.
        if dx(idx) == 0 && dy(idx) == 0
    
            pX = x(1,idx) + radius(idx)*cos(tCir);
            pY = y(1,idx) + radius(idx)*sin(tCir);
            cPS = polyshape(pX,pY);
            minZ = min(z(:,idx));
            maxZ = max(z(:,idx));
        else
        
            %find orthogonal vectors @ endpoints
            aVp1 = [aV(idx,2),-aV(idx,1)];
            aVp2 = [-aV(idx,2),aV(idx,1)];
            bVp1 = [bV(idx,2),-bV(idx,1)];
            bVp2 = [-bV(idx,2),bV(idx,1)];
        
            aVp1 = aVp1./norm(aVp1);
            aVp2 = aVp2./norm(aVp2);
            bVp1 = bVp1./norm(bVp1);
            bVp2 = bVp2./norm(bVp2);
        
            %from each endpoint, use radius to find vertices of the rectangle
            x1 = x(1,idx) + radius(idx)*aVp1(1);
            y1 = y(1,idx) + radius(idx)*aVp1(2);
            x2 = x(1,idx) + radius(idx)*aVp2(1);
            y2 = y(1,idx) + radius(idx)*aVp2(2);  
            x3 = x(2,idx) + radius(idx)*bVp1(1);
            y3 = y(2,idx) + radius(idx)*bVp1(2);
            x4 = x(2,idx) + radius(idx)*bVp2(1);
            y4 = y(2,idx) + radius(idx)*bVp2(2);
    
    
            %calculate set of orthgonal vectors
            ZOrtho = oVz(aV(idx,:),XOrtho,YOrtho);
        
            %unit-ify the orthgonal vectors
            uov = [XOrtho,YOrtho,ZOrtho]./sqrt(XOrtho.^2 + YOrtho.^2 + ZOrtho.^2);
    
              
            %donot re unit-fy, you only want the horizontal component, not the
            %renormalized horizontal component
        
            %using only the X and Y components, find circle coods in plane of
            %interest
            xaC = x(1,idx) + uov(:,1).*radius(idx);
            yaC = y(1,idx) + uov(:,2).*radius(idx);
            zaC = z(1,idx) + uov(:,3).*radius(idx);
    
            xbC = x(2,idx) + uov(:,1).*radius(idx);
            ybC = y(2,idx) + uov(:,2).*radius(idx);
            zbC = z(2,idx) + uov(:,3).*radius(idx);
    
            minZ(idx) = min([zaC; zbC]);
            maxZ(idx) = max([zaC; zbC]);
        
    %         plot(x1,y1,'k.')
    %         hold on
    %         plot(x2,y2,'ko')
    %         plot(x3,y3,'k^')
    %         plot(x4,y4,'k*')
    %     
    %         plot(x(1,idx)',y(1,idx)','ro')
    %         plot(x(2,idx)',y(2,idx)','bo')
    %         
    %         plot3(x(:,idx),y(:,idx),z(:,idx))
    %         plot(xaC,yaC,'m.')
    %         plot(xbC,ybC,'g.')
    %         axis equal
            
            %assymble total package
            rX = [x1; x2; x3; x4];
            rY = [y1; y2; y3; y4];
    
            %test for circle parts in polygon
    
       
    
            partsPS = [polyshape(xaC,yaC),...
                polyshape(rX,rY),...
                polyshape(xbC,ybC)];
    
            cPS = union(partsPS);
            cPS = reducepoly(cPS.Vertices);
            
       
        end
        
        pSV(idx) = cPS;%save polygon
    
        if rand < 0.04
            fprintf('completed %i \n',round((idx/noCylinders)*100))
        end
    
    end
    T.pSV = pSV;
    T.minZ = minZ;
    T.maxZ = maxZ;
end