function [dist, dtwM, path] = OBEDTW(q,r,t)
% This function implements the OpenBeginEndDTW with no warping constraints
% For details, please check the paper "Prefix and Suffix Invariant 
% Dynamic Time Warping"
% --------------------------------------------------
% Input parameters:
% - q and r - (row) vectors of real numbers describing the time series
% - t is the drop off ratio at the start and the end of time series signal
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

% q = [1,3,4,3,4];
% r=[0,2,2,2,4];
% t = 0.2;
    n = length(q);
    m = length(r);

% DTW matrix initialization
    
    dtwM = inf(n+1,m+1);
 
    dtwM(1:round(t*n)+1,1:round(t*m)+1) = 0;
    
    dtwM(n+1-round(n*t):n+1, m+1-round(m*t):m+1) = 0;
 
    path = [];
    
% Loops related to the DTW's recurrence relation
% calculate dtw_prefix
     %calculate the first part of matrix
     for i = 2:round(n*t) +1
         for j = round(t*m)+2 : m+1
             
             dtwM(i,j) = (q(i-1)-r(j-1))^2 + ...
                 min(dtwM(i-1,j-1), min(dtwM(i-1,j), dtwM(i,j-1)));
             
         end
         
     end
     %calculate the second part of matrix
     for i = round(t*n)+2 : n+1-round(t*n)
         for j = 2 : m+1
             
             dtwM(i,j) = (q(i-1)-r(j-1))^2 + ...
                 min(dtwM(i-1,j-1), min(dtwM(i-1,j), dtwM(i,j-1)));
             
         end
         
     end
     %calculate the third part of matrix
     for i = n+2-round(t*n) : n+1
         for j = 2 : m+1-round(m*t)
             
             dtwM(i,j) = (q(i-1)-r(j-1))^2 + ...
                 min(dtwM(i-1,j-1), min(dtwM(i-1,j), dtwM(i,j-1)));
             
         end
         
     end   
    
    
    % Final distance

    
     bottomright_line = min(dtwM((n+1-round(n*t)):n+1,m+1-round(m*t)));
     rightbottom_line = min(dtwM((n+1-round(n*t)),m-round(m*t):m+1));
     
     dist = min(bottomright_line,rightbottom_line);

    % in the case we are interested in the optimal alignment
    if (nargout > 2)
        path = getPath(dtwM,t);
    end
    csvwrite('OBEDTW_matrix.csv',dtwM);
    [x,y] = size(dtwM);
 
    M1 = dtwM(x-round((x-1)*t):x,y-round((y-1)*t));
    M2 = dtwM(x-round((x-1)*t),y-round((y-1)*t):y);
    temp_col = min(M1,[],1);
    temp_row = min(M2,[],2);
    smallest = min(temp_col,temp_row);
   
    [r,c] = find(dtwM==smallest);
    xixi = r(end);
    haha = c(end);
%     disp(xixi);
%     disp(haha);
   

    
end

% This function traverses the DTW matrix in order to recover the alignment
% obtained by the algorithm. In our paper, we used it to plot figures

function path = getPath(dtwM,t)
     [x,y] = size(dtwM);
     M1 = dtwM(x-round((x-1)*t):x,y-round((y-1)*t));
     M2 = dtwM(x-round((x-1)*t),y-round((y-1)*t):y);
     temp_col = min(M1,[],1);
     temp_row = min(M2,[],2);
     
     smallest = min(temp_col,temp_row);
     
     [r,c] = find(dtwM==smallest);

     m = r(end);
     n = c(end);
%      disp(m);
%      disp(n);
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
     
    
     
 
 end

