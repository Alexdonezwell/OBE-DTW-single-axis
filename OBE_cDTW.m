function [dist, dtwM, path] = OBE_cDTW(q,r,warpingRate,t)
% This function implements the OpenBeginEndDTW with no warping constraints
% For details, please check the paper "Prefix and Suffix Invariant 
% Dynamic Time Warping"
% --------------------------------------------------
% Input parameters:
% - q and r - (row) vectors of real numbers describing the time series
% - t is the drop off ratio at the start and the end of time series signal
% - warpingRate - the relative length of the sakoe-chiba warping window
% --------------------------------------------------
% Output parameters:
% - dist - the final distance estimate
% - dtwM - the matrix (of cumulative cost) used to calculate the distance
% - path - a vector representation of the alignment path
% --------------------------------------------------
% Examples:
% test case: q=[1,3,4,3,4]    r=[0,2,2,2,4]
% [dist, dtwM, path] = OBEDTW(q, r,.2)
% dtwM =

%      0     0   Inf   Inf   Inf   Inf
%      0     0     1     2     3    12
%    Inf     9     1     2     3     4
%    Inf    25     5     5     6     3
%    Inf    34     6     6     6     4
%    Inf    50    10    10    10     0


% [dist, dtwM, path] = OBEDTW(q, r,.4)
% dtwM =
 
%      0     0     0   Inf   Inf   Inf
%      0     0     0     1     2    11
%      0     0     0     1     2     3
%    Inf    16     4     4     5     2
%    Inf    25     5     5     0     0
%    Inf    41     9     9     0     0
% --------------------------------------------------


    n = length(q);
    m = length(r);

% DTW matrix initialization
    
    dtwM = inf(n+1,m+1);
    dtwM(1:round(t*n)+1,1:round(t*m)+1) = 0;
    dtwM(n+1-round(n*t):n+1, m+1-round(m*t):m+1) = 0;
    
% defines the number of cells in the warping window
    winLen = ceil(n * warpingRate);
    W_W = max(winLen, abs(n-m));

% 
    path = [];
    
% Loops related to the DTW's recurrence relation
% calculate dtw_prefix
     %calculate the first part of matrix
     for i = 2:round(n*t) +1
         beg_fake1 = max(2, i-W_W);
         end_fake1 = max(round(m*t) +1, i+W_W);
         start1 = max(beg_fake1,round(t*m)+2);
         end1 = min(end_fake1,m+1);
         for j = start1:end1
             
             dtwM(i,j) = (q(i-1)-r(j-1))^2 + ...
                 min(dtwM(i-1,j-1), min(dtwM(i-1,j), dtwM(i,j-1)));
             
         end
         
     end
     %calculate the second part of matrix
     for i = round(t*n)+2 : n+1-round(t*n)
         end_fake2 = min(n+1-round(t*m), i+W_W);
         start2 = max(i-W_W,2);
         end2 = min(end_fake2,m+1);
         for j = start2:end2
             
             dtwM(i,j) = (q(i-1)-r(j-1))^2 + ...
                 min(dtwM(i-1,j-1), min(dtwM(i-1,j), dtwM(i,j-1)));
             
         end
         
     end
     %calculate the third part of matrix
     for i = n+2-round(t*n) : n+1
         beg_fake3 = max(n+2-round(t*m), i-W_W);
         end_fake3 = min(m+1, i+W_W);
         start3 = max(beg_fake3,2);
         end3 = min(end_fake3,m+1-round(m*t));
         for j = start3:end3
             
             dtwM(i,j) = (q(i-1)-r(j-1))^2 + ...
                 min(dtwM(i-1,j-1), min(dtwM(i-1,j), dtwM(i,j-1)));
             
         end
         
     end
     %Deal with the Boundary conditions
%      for i = n+1-round(t*n) 
%          for j = m+1-round(m*t):m+1
%              dtwM(i,j) = (q(i-1)-r(j-1))^2 + ...
%                  min(dtwM(i-1,j-1), min(dtwM(i-1,j), dtwM(i,j-1)));
%              
%          end
%          
%      end
%      
%      for i = n+1-round(t*n):n+1 
%          for j = m+1-round(m*t)
%              dtwM(i,j) = (q(i-1)-r(j-1))^2 + ...
%                  min(dtwM(i-1,j-1), min(dtwM(i-1,j), dtwM(i,j-1)));
%              
%          end
%          
%      end
    
    
    % Final distance

    
     bottomright_line = min(dtwM((n+1-round(n*t)):n+1,m+1-round(m*t)));
     rightbottom_line = min(dtwM((n+1-round(n*t)),m-round(m*t):m+1));
     
     dist = min(bottomright_line,rightbottom_line);

    % in the case we are interested in the optimal alignment
    if (nargout > 2)
        path = getPath(dtwM);
    end
    csvwrite('OBE_cDTW_matrix.csv',dtwM);
end

% This function traverses the DTW matrix in order to recover the alignment
% obtained by the algorithm. In our paper, we used it to plot figures
function path = getPath(dtwM)

    [~,idxMinRow] = min(dtwM(end,:));
        
    n = size(dtwM,2);
    m = idxMinRow;

    path = [];
    
    while (n > 1) %first column or row reached
        
        path = [path; n, m];
        
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

    end
    
    path = path - 1;

end

