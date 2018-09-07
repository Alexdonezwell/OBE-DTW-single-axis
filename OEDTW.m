function [dist, dtwM, path] = OEDTW(q, r)
% This function implements the OpenEndDTW with no warping constraints
% For details, please check the paper "Prefix and Suffix Invariant 
% Dynamic Time Warping"

% Input parameters:
% - q and r - (row) vectors of real numbers describing the time series

% Output parameters:
% - dist - the final distance estimate
% - dtwM - the matrix (of cumulative cost) used to calculate the distance
% - path - a vector representation of the alignment path

    n = length(q);
    m = length(r);

    % DTW matrix initialization
    dtwM = inf(n+1,m+1);
    dtwM(1,1) = 0;
    
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
    
    % Final distance
    dist = min(dtwM(end,:));

    % in the case we are interested in the optimal alignment
    if (nargout > 2)
        path = getPath(dtwM);
    end
    
end

% This function traverses the DTW matrix in order to recover the alignment
% obtained by the algorithm. In our paper, we used it to plot figures
function path = getPath(dtwM)

    [~,idxMinRow] = min(dtwM(end,:));
        
    n = size(dtwM,2);
    m = idxMinRow;

    path = [];
    
    while (m > 1 || n > 1)
        
        path = [path; n, m];
               
        if (m > 1 && n > 1)
            
            [~, minPos] = min([dtwM(n-1,m-1), ...
                dtwM(n,m-1), dtwM(n-1,m)]);

            if (minPos == 1)
                m = m-1;
                n = n-1;
            else
                if (minPos == 2)
                    m = m-1;
                else
                    n = n-1;
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

