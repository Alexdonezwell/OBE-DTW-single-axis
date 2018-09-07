function [dist, dtwM, path] = cPsiDTW(q, r, percToOpen, warpingRate, bsf)
% This function implements the psi-DTW with the warping constraint
% Note: this function considers that both time series are equal length
% For details, please check the paper "Prefix and Suffix Invariant 
% Dynamic Time Warping"

% Input parameters:
% - q and r - (row) vectors of real numbers describing the time series
% - percToOpen - the relative relaxing factor, i. e., a [0 1] number
% - warpingRate - the relative length of the sakoe-chiba warping window
% - bsf - the best-so-far used by the 1NN algorithm (optional)

% Output parameters:
% - dist - the final distance estimate
% - dtwM - the matrix (of cumulative cost) used to calculate the distance
% - path - a vector representation of the alignment path

    % The time series have the same length
    n = length(q);
    
    % defines the number of cells in the warping window
    winLen = ceil(n * warpingRate);
    % defines the number of cells to "relax"
    tolerance = ceil(n * percToOpen);

    % DTW matrix initialization
    dtwM = inf(n+1,n+1);
    dtwM(1:tolerance+1,1) = 0;
    dtwM(1,1:tolerance+1) = 0;
    
    path = [];
    
    % Loops related to the DTW's recurrence relation
    for i = 2 : n+1
        
        begObs = max(2, i-winLen);
        endObs = min(n+1, i+winLen);
        
        for j = begObs : endObs
            
            dtwM(i,j) = (q(i-1)-r(j-1))^2 + ...
                min(dtwM(i-1,j-1), min(dtwM(i-1,j), dtwM(i,j-1)));
            
        end
        
        % Early abandon
        if (nargin == 5 && min(dtwM(i,:)) > bsf)
            dist = inf;
            return;
        end
        
    end
    
    % finds the min value in the "relaxed" region
    minLine = min(dtwM(end-tolerance:end,end));
    minColumn = min(dtwM(end,end-tolerance:end));
    dist = min(minLine, minColumn);
    
    % in the case we are interested in the optimal alignment
    if (nargout > 2)
        path = getPath(dtwM, tolerance, tolerance);
    end
    csvwrite('PSIcDTW_matrix.csv',dtwM);
end

% This function traverses the DTW matrix in order to recover the alignment
% obtained by the algorithm. In our paper, we used it to plot figures
function path = getPath(dtwM, toleranceN, toleranceM)

    [vRow,idxMinRow] = min(dtwM(end-toleranceN:end,end));
    [vColumn,idxMinColumn] = min(dtwM(end,end-toleranceM:end));

    if(vRow < vColumn)
        
        m = size(dtwM,1) - toleranceM + idxMinRow - 1;
        n = size(dtwM,2);
        
    else
        
        n = size(dtwM,2) - toleranceN + idxMinColumn - 1;
        m = size(dtwM,1);
        
    end

    path = [];
    
    while (m > 1 && n > 1)
        
        path = [path; m,n];
        
        [~, minPos] = min([dtwM(m-1,n-1), ...
            dtwM(m,n-1), dtwM(m-1,n)]);

        if (minPos == 1)
            m = m-1;
            n = n-1;
        else
            if (minPos == 2)
                n = n-1;
            else
                m = m-1;
            end
        end

    end
    
    while (m == 1 && n > toleranceN + 1)        
        path = [path; m,n];
        n = n-1;        
    end
    
    while (n == 1 && m > toleranceM + 1)        
        path = [path; m,n];
        m = m-1;        
    end
    
    path = path - 1;

end

