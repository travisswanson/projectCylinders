function T =  readFiles(FN,FD,T)
    cd(FD)
    for idx = 1:length(FN)
        fprintf(['reading...' FN{idx} '\n'])
        [DATA, ~, ~] = tblread(FN{idx},',');
    
    
        x = [DATA(:,3) DATA(:,6)]';
        y = [DATA(:,4) DATA(:,7)]';
        z = [DATA(:,5) DATA(:,8)]';
        radius = DATA(:,9);
        
        noCylinders = size(DATA,1);
        cLength = DATA(:,12);
        
        dx = DATA(:,6)- DATA(:,3);
        dy = DATA(:,7)- DATA(:,4);
        dz = DATA(:,8)- DATA(:,5);
        
        theta = atand(dz./sqrt(dx.^2 + dy.^2));
        
        BO = DATA(:,20);
        maxBO = max(BO);
        
        T(idx).x = x;
        T(idx).y = y;
        T(idx).z = z;
        T(idx).radius = radius;
        T(idx).noCylinders = noCylinders;
        T(idx).cLength = cLength;
        T(idx).dx = dx;
        T(idx).dy = dy;
        T(idx).dz = dz;
        T(idx).theta = theta;
        T(idx).BO = BO;
        T(idx).maxBO = maxBO;

    end
end