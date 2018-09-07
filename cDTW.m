function [dist, dtwM, path] = cDTW(q, r, warpingRate, bsf)
% This function implements the DTW with the warping constraint
% Note: this function considers that both time series are equal length

% Input parameters:
% - q and r - (row) vectors of real numbers describing the time series
% - warpingRate - the relative length of the sakoe-chiba warping window
% - bsf - the best-so-far used by the 1NN algorithm (optional)

% Output parameters:
% - dist - the final distance estimate
% - dtwM - the matrix (of cumulative cost) used to calculate the distance
% - path - a vector representation of the alignment path

    % The time series have the same length
    n = length(q);

    % DTW matrix initialization
    dtwM = inf(n+1,n+1);
    dtwM(1,1) = 0;

    % defines the number of cells in the warping window
    winLen = ceil(n * warpingRate);
    
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
        if (nargin == 4 && min(dtwM(i,:)) > bsf)
            dist = inf;
            return;
        end
        
    end
    
    % Final distance value
    dist = dtwM(end,end);

    % in the case we are interested in the optimal alignment
    if (nargout > 2)
        path = getPath(dtwM);    
    end
    csvwrite('cDTW_matrix.csv',dtwM);
end

% This function traverses the DTW matrix in order to recover the alignment
% obtained by the algorithm. In our paper, we used it to plot figures
function path = getPath(dtwM)

    [m,n] = size(dtwM);
    path = [];
    
    while (m > 1 || n > 1)
        
        path = [path; m,n];
               
        if (m > 1 && n > 1)
            
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
            
        else
            
            if (m > 1) % n == 1                
                m = m-1;                
            else
                n = n-1;
            end
            
        end
                
    end
    
    path = path - 1;

end

