function [dist, dtwM, path] = psiDTW(q, r, percToOpen, bsf)
% This function implements the psi-DTW with no warping constraints
% For details, please check the paper "Prefix and Suffix Invariant 
% Dynamic Time Warping"

% Input parameters:
% - q and r - (row) vectors of real numbers describing the time series
% - percToOpen - the relative relaxing factor, i. e., a [0 1] number
% - bsf - the best-so-far used by the 1NN algorithm (optional)

% Output parameters:
% - dist - the final distance estimate
% - dtwM - the matrix (of cumulative cost) used to calculate the distance
% - path - a vector representation of the alignment path

    n = length(q);
    m = length(r);
    
    % defines the number of cells to "relax"
    toleranceN = ceil(n * percToOpen);
    toleranceM = ceil(m * percToOpen);

    % DTW matrix initialization
    dtwM = inf(n+1,m+1);
    dtwM(1:toleranceN+1,1) = 0;
    dtwM(1,1:toleranceM+1) = 0;
    
    path = [];
    
    % Loops related to the DTW's recurrence relation
    for i = 2 : n+1
        for j = 2 : m+1
            
            dtwM(i,j) = (q(i-1)-r(j-1))^2 + ...
                min(dtwM(i-1,j-1), min(dtwM(i-1,j), dtwM(i,j-1)));
            
        end
        
        % Early abandon
        if (nargin == 4 && min(dtwM(i,:)) > bsf)
            dist = inf;
            return;
        end
        
    end
    
    % finds the min value in the "relaxed" region
    minRow = min(dtwM(end-toleranceN:end,end));
    minColumn = min(dtwM(end,end-toleranceM:end));
    dist = min(minRow, minColumn);

    % in the case we are interested in the optimal alignment
    if (nargout > 2)
        path = getPath(dtwM, toleranceN, toleranceM);
    end
    csvwrite('PSIDTW_matrix.csv',dtwM);
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

