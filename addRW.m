function newData = addRW(data, RWLen, seed)

    if (nargin < 3)
        rng(1);
    else
        rng(seed)
    end

    newData = zeros(size(data,1), size(data,2) + ceil(RWLen/2)*2);
    
    for i = 1 : size(data,1)
        
        minV = min(data(i,:));
        maxV = max(data(i,:));
        
        preTS = newRW(ceil(RWLen/2), minV, maxV, data(i,1));
        % Fix to match the end of RW to the beginning of the original
        preTS = fliplr(preTS);
                
        posTS = newRW(ceil(RWLen/2), minV, maxV, data(i,end));
        
        newData(i,:) = [ preTS, data(i,:), posTS];
        
    end

end

function rw = newRW(len, minV, maxV, referenceValue)

    rw = rand(1,len) - 0.5;
    %newRW = cumsum(rw);
    
    rw(1) = rw(1) + referenceValue;    
    for i = 2 : len
        
        nextVal = rw(i) + rw(i-1);
        
        % if it cross the limits
        if (nextVal > maxV || nextVal < minV)
            nextVal = rw(i) - rw(i-1);
        end
        
        rw(i) = nextVal;
        
    end

end