function T = cylinderUnion(T,chunkSize,itCount)

    minZ = T.minZ;
    pSV = T.pSV;
    noCylinders = T.noCylinders;

    [sminZ,sIdx]=sort(minZ,'descend');
    spSV = pSV(sIdx);
    
    kdx = 1;

    uPS = spSV(1);
    rawArea = spSV(1).area;

    chunkRawArea = zeros(1,ceil(noCylinders/chunkSize));
    chunkMinZ = chunkRawArea;
    lMinZ = 1;

    for idx = 2:noCylinders 
        rawArea = rawArea + spSV(idx).area;
        uPS = union(uPS,spSV(idx));
        
    
        if mod(idx,chunkSize) == 0 && ~(idx == noCylinders)
            chunkUnion(kdx) = uPS;
            uPS = spSV(idx+1);
           
            chunkRawArea(kdx) = rawArea;
            rawArea = 0;
            chunkMinZ(kdx) = mean(sminZ(lMinZ:idx));
            lMinZ = idx;
            idx = idx + 2;
            kdx = kdx + 1;
        elseif idx == noCylinders
       
            chunkRawArea(kdx) = rawArea;
            rawArea = 0;
            chunkUnion(kdx) = uPS;
            chunkMinZ(kdx) = mean(sminZ(lMinZ:idx));
        end
        if rand < 0.05
            fprintf('Tree %i - Chunk union - just completed %i \n',itCount, idx)
        end
    end
    
    
    
    sF = zeros(1,length(chunkUnion));
    
    ultimateUnion = chunkUnion(1);
    sF(1) = chunkUnion(1).area/chunkRawArea(1);
    
    for idx = 2:length(chunkUnion)
        ultimateUnion = union(ultimateUnion,chunkUnion(idx));
        sF(idx) = ultimateUnion.area/sum([chunkRawArea(1:idx)]);
        fprintf('Tree %i - Ultiamte union - just completed %i \n',itCount, idx)
    end

    T.ultimateUnion = ultimateUnion;
    T.shadowFraction = sF;
    T.chunkUnion = chunkUnion;
    T.chunkRawArea = chunkRawArea;
    T.chunkMinZ = chunkMinZ;

end