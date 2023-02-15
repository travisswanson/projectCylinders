clear variables
T = initializeStructure();

FD = '/home/travis/data/';

FN = dir('/home/travis/data/*.csv');
FN = {FN.name};

T =  readFiles(FN,FD,T);
cd('/home/travis/')

parpool(20)

% Project Cylinders
parfor idx = 1:length(T)
    T(idx) = projectCylinders(T(idx));
end


% Union Cylinders
chunkSize = 200;

parfor idx = 1:length(T)
    T(idx) = cylinderUnion(T(idx),chunkSize,idx); 
end


%% VISUALIZE 
for idx = 1:length(T)
    subplot(4,4,idx)
        ha = plot(T(idx).ultimateUnion);
        title(FN(idx))
        xlim([min(T(idx).x(:)) max(T(idx).x(:))])
        ylim([min(T(idx).y(:)) max(T(idx).y(:))])
         hold on
        axis equal
        axis off
      
       
end
set(gcf,'color','w')

%% shadow fraction
for idx = 1:length(T)
    subplot(4,4,idx)
        ha = plot(T(idx).ultimateUnion);
        xlim([min(T(idx).x(:)) max(T(idx).x(:))])
        ylim([min(T(idx).y(:)) max(T(idx).y(:))])
         hold on
        axis equal
        axis off
       
end
set(gcf,'color','w')

cm = zeros(13,3);
for idx = 1:length(T)
    subplot(4,4,idx)
    plot(cumsum(T(idx).chunkRawArea.*T(idx).shadowFraction),T(idx).chunkMinZ,'color',cm(idx,:),'linewidth',1)
    hold on
    plot(cumsum(T(idx).chunkRawArea.*(1-T(idx).shadowFraction)),T(idx).chunkMinZ,'--','color',cm(idx,:),'linewidth',1)
    xlabel('Canopy area m^{2}')
    ylabel('Elevation m')
    axis tight
end
legend({'Exposed Canopy','Shadowed Canopy'})

%% Make a PDF

parfor idx = 1:length(T)
        f = figure('visible', 'off');
        ha = plot(T(idx).ultimateUnion);
        title(FN(idx))
        xlim([min(T(idx).x(:)) max(T(idx).x(:))])
        ylim([min(T(idx).y(:)) max(T(idx).y(:))])
        set(gcf,'color','w')
        hold on
        axis equal
        axis off
        print([num2str(idx) '_' FN{idx} '.png'],'-dpng')
        close(f)
 
end





%%

% x = x(:,1);
% y = y(:,1);
% z = z(:,1);
% 
% x = [x(1:end) x(2:)];
for idx = 1:10:noCylinders
    plot3(x(:,idx),y(:,idx),z(:,idx),'k','linewidth',DATA(idx,9)*50)
    hold on
end

%% Branch angle 



histogram(theta)
xlabel('Branch Dip Angle (°)')

%% Branch dip angle as a function of elevation

cm = viridis(maxBO);


for idx = 1:maxBO
    li = BO == idx;
    [F,X] = ecdf(theta(li));
    
    plot(X,F,"Color",cm(idx,:),'LineWidth',2)
    hold on
    ltxt{idx} = ['Branch Order ' num2str(idx)];
end
ylabel('Cum. Prob.')
xlabel('Branch Dip Angle (°)')
legend(ltxt)
set(gcf,'color','w')
set(gca,'fontsize',14)

%%



%% Intersect

iPS = pSV(1);
kdx =1;
for idx = 1:noCylinders
    
    iPS = xor(iPS,pSV(idx));

    fprintf('just completed %i \n',idx)
    if mod(idx,500) == 0
        chunkIntersect(kdx) = iPS;
        iPS = pSV(idx+1);
        idx = idx + 2;
        kdx = kdx + 1;
    
    end
end

%%
ultimateIntersect = chunkIntersect(1);
for idx = 2:length(chunkIntersect)
    ultimateIntersect = xor(ultimateIntersect,chunkIntersect(idx));
    fprintf('just completed %i \n',idx)
end

%% Find how much of a cylinder is unique?

uFA = nan(1,noCylinders);
nuFA = nan(1,noCylinders);

for idx = 1%:noCylinders

    uF = intersect(ultimateIntersect,pSV(idx));
    

end

% plot(uPS);
% axis equal

%% Better way

%Calculate centroids, then pair-wise distance, then find all centroids
%within range of max cylinder size.

%find overlapping bounding boxes (very cautious!)
xb = nan(noCylinders,2);
yb = xb;
for idx = 1:noCylinders
    [xb(idx,:),yb(idx,:)]=pSV(idx).boundingbox;
end

xb = [xb(:,1) xb(:,2) xb(:,2) xb(:,1) xb(:,1)];
yb = [yb(:,1) yb(:,1) yb(:,2) yb(:,2) yb(:,1)];

%% The way it should have been done all along


%%
plot(sF.*cumsum(chunkRawArea),chunkMinZ,'k','linewidth',1)
hold on
plot((1-sF).*cumsum(chunkRawArea),chunkMinZ,'k--','linewidth',1)
ylabel('Elevation m')
xlabel('Projected Canopy Area m^{2}')
ylim([0 15])
set(gca,'fontsize',14)
set(gcf,'color','w')
legend({'Exposed Canopy','Shadowed Fraction'})
%%
tic


sFV = nan(noCylinders,1);

parpool(40)

parfor kdx = 1:noCylinders
    idxList = 1:noCylinders;

    idx2test = kdx;
    hitList = false(noCylinders,1);

    for idx = 1:noCylinders

        if idx == idx2test
            continue
        end
        tmp = pSV(idx).Vertices;
        [in, on] = inpolygon(xb(idx2test,:),yb(idx2test,:),tmp(:,1),tmp(:,2));
        hitList(idx) = any(in|on);

    end
    
    idx2Consider = idxList(hitList);

    %find all cylinders with end points above minimum elevation of COI
    %liElev = (min(z(:,idx2test))) <= (min(z(:,idx2Consider))');
    
    %lowest possible elevation on branch of interest
    elBOI= minZ(idx2test);% -radius(idx2test);
    elHitLow = minZ(idx2Consider);% - radius(idx2Consider)';
    elHitHigh = maxZ(idx2Consider);%+ radius(idx2Consider)';

    liElev = (elBOI <= elHitLow) | (elBOI <= elHitHigh);

    %only these cylinders can shade the COI
    if nnz(liElev) == 0
        sF = 0;
    else
        chunk = union(pSV(idx2Consider(liElev)));
        shade = intersect(chunk,pSV(idx2test));
        sF = shade.area/pSV(idx2test).area;

    end

    sFV(kdx) = sF;


    
% smpIdx = [idx2test idx2Consider];
% for idx = 1:length(smpIdx)
% 
%         if smpIdx(idx) ~= idx2test
%             plotCylinder(x(:,smpIdx(idx)),y(:,smpIdx(idx)),z(:,smpIdx(idx))...
%                 ,dx(smpIdx(idx)),dy(smpIdx(idx)),dz(smpIdx(idx)),radius(smpIdx(idx)),'g')
%             hold on
%         else
%             plotCylinder(x(:,smpIdx(idx)),y(:,smpIdx(idx)),z(:,smpIdx(idx))...
%             ,dx(smpIdx(idx)),dy(smpIdx(idx)),dz(smpIdx(idx)),radius(smpIdx(idx)),'r')
%             hold on
%         end
% 
%  
% end


title(['Shadowed Fraction ' num2str(sF) ])
    
    if rand < 0.01
        fprintf('just completed %i \n',kdx)
    end
end
toc


 
%% Visualize area as a function of elevation
mZ = mean(z,1);
rawA = [pSV.area];
[~,sIdx]=sort(mZ,'descend');
mZ = mZ(sIdx);
rawA = rawA(sIdx);
sF = sFV(sIdx);
plot(cumsum(rawA),mZ,'k-','linewidth',1)
hold on
plot(cumsum(rawA(:).*(1-sF(:))),mZ,'r-','linewidth',1)
box on 
xlabel('Area m^{2}')
ylabel('Elevation m')
set(gca,'fontsize',14)
set(gcf,'color','w')
legend({'Total Cum. Branch Area','Non-shadowed Cum. Branch Area'})
%% brute force way

idxList = 1:noCylinders;
idxSplit = ceil(linspace(1,noCylinders,500));

idx2save = idxSplit(idx-1):idxSplit(idx);
ppSV = pSV;
ppSV(idx2save) = [];

uPS = ppSV(1);
kdx =1;
for idx = 1:length(ppSV)
    
    uPS = union(uPS,ppSV(idx));
    fprintf('just completed %i \n',idx)
    if mod(idx,500) == 0
        chunkUnion(kdx) = uPS;
        uPS = pSV(idx+1);
        idx = idx + 2;
        kdx = kdx + 1;
    
    end
end


puUnion = chunkUnion(1);
for idx = 2:length(chunkUnion)
    puUnion = union(puUnion,chunkUnion(idx));
    fprintf('just completed %i \n',idx)
end



for idx = 1:len
end

tic


sFV = nan(noCylinders,1);

% parpool(40)

parfor kdx = 1:noCylinders
    

    idx2test = kdx;
    hitList = false(noCylinders,1);

    for idx = 1:noCylinders

        if idx == idx2test
            continue
        end

        [in, on] = inpolygon(xb(idx2test,:),yb(idx2test,:),xb(idx,:),yb(idx,:));
        hitList(idx) = any(in|on);

    end
    
    idx2Consider = idxList(hitList);

    %find all cylinders with end points above minimum elevation of COI
    liElev = (min(z(:,idx2test))) <= (min(z(:,idx2Consider))');
    liElev = (min(z(:,idx2test))-radius(idx2test)) <= (min(z(:,idx2Consider)) + radius(idx2Consider)');

    %only these cylinders can shade the COI
    if nnz(liElev) == 0
        sF = 0;
    else
        chunk = union(pSV(idx2Consider(liElev)));
        shade = intersect(chunk,pSV(idx2test));
        sF = shade.area/pSV(idx2test).area;

    end

    sFV(kdx) = sF;


    
% smpIdx = [idx2test idx2Consider];
% for idx = 1:length(smpIdx)
% 
%         if smpIdx(idx) ~= idx2test
%             plotCylinder(x(:,smpIdx(idx)),y(:,smpIdx(idx)),z(:,smpIdx(idx))...
%                 ,dx(smpIdx(idx)),dy(smpIdx(idx)),dz(smpIdx(idx)),radius(smpIdx(idx)),'g')
%             hold on
%         else
%             plotCylinder(x(:,smpIdx(idx)),y(:,smpIdx(idx)),z(:,smpIdx(idx))...
%             ,dx(smpIdx(idx)),dy(smpIdx(idx)),dz(smpIdx(idx)),radius(smpIdx(idx)),'r')
%             hold on
%         end
% 
%  
% end


title(['Shadowed Fraction ' num2str(sF) ])
    
    if rand < 0.01
        fprintf('just completed %i \n',kdx)
    end
end
toc
 
%% SCRAPS


% if idx == 1
%     pU = union(pSV(smpIdx(1)),pSV(smpIdx(2)));
% end
%  plot(pU)
% 
% hold on
% plot(pSV(idx2test))
% % axis equal
% title(['Shadowed Fraction ' num2str(sF) ])
% 
% for idx = 1:length(smpIdx)
%         if smpIdx(idx) ~= idx2test
%     plotCylinder(x(:,smpIdx(idx)),y(:,smpIdx(idx)),z(:,smpIdx(idx))...
%         ,dx(smpIdx(idx)),dy(smpIdx(idx)),dz(smpIdx(idx)),radius(smpIdx(idx)),'g')
%     hold on
%         else
%                 plotCylinder(x(:,smpIdx(idx)),y(:,smpIdx(idx)),z(:,smpIdx(idx))...
%         ,dx(smpIdx(idx)),dy(smpIdx(idx)),dz(smpIdx(idx)),radius(smpIdx(idx)),'r')
%         end
% 
%  
% end

% 
%     if length(smpIdx) == 1
%         %in the case of only one hit cylinder, it's totally exposed
%         sF = 0;
%     else
%         for idx = 1:length(sIdx)
%     
%             if idx == 1
%     
%                 if smpIdx(1) == idx2test 
%                     iB = intersect(pSV(idx2test),pSV(smpIdx(2)));
%                     sF = (iB.area)./pSV(idx2test).area;
%                     break
%                 elseif smpIdx(2) == idx2test
%                     iB = intersect(pSV(idx2test),pSV(smpIdx(1)));
%                     sF = (iB.area)./pSV(idx2test).area;
%                     break
%                 else
%                     pU = union(pSV(smpIdx(1)),pSV(smpIdx(2)));
%                 end
%     
%             else
%     
%                 if smpIdx(idx) == idx2test
%                     iB = intersect(pSV(idx2test),pU);
%                     sF = iB.area./pSV(idx2test).area;
%                     break
%                 else
%                     pU = union(pU,pSV(smpIdx(idx)));
%                 end
%     
%             end
%        
%         end
