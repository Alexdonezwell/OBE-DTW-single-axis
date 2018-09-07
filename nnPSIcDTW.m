function [accuracy, predictions, nearestEx, distMatrix] = ...
            nnPSIcDTW(data, percToRelax)
                                   
    nTr = length(data.train);
    nTe = length(data.test);
    
    predictions = zeros(nTe,1);
    nearestEx = zeros(nTe,1);
    
    distMatrix = zeros(nTe, nTr);
    accuracy = 0;
    
    if (percToRelax <= 0.1)
        thisWLen = 0.1;
    else
        thisWLen = percToRelax;
    end

    %z-norm test
    for i = 1 : nTe
        data.test(i).data = ...
            znorm(data.test(i).data);
    end
    
    %z-norm train
    for i = 1 : nTr
        data.train(i).data = znorm(data.train(i).data);
    end
    
    for i = 1 : nTe
        
        display([datestr(now), ' - Te #', num2str(i)])
        
        best = inf;
        
        for j = 1 : nTr

			% Resize the training subsequences to the length of the query
            thisResizedTrain = imresize(data.train(j).data, ...
                size(data.test(i).data));
            
            thisDist = cPsiDTW(...
                data.test(i).data', ...
                thisResizedTrain', percToRelax, thisWLen, best);
            
            distMatrix(i,j) = thisDist;
            
            if (thisDist < best)
                best = thisDist;
                predictions(i) = data.train(j).label;
                nearestEx(i) = j;
            end
            
        end
        
        if (predictions(i) == data.test(i).label)
            accuracy = accuracy + 1;
        end

    end
    
    accuracy = accuracy / nTe;
    csvwrite('PSIcDTW_matrix.csv',distMatrix);
    
end

