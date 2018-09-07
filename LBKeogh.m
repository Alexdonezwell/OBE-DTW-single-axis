function [LB] = LBKeogh(q,r,relativeW)

    len = length(q);
    
    w = ceil(len * relativeW);
    LB = 0;
    
    for i = 1 : len

        thisSubSeq = q(max(1,i-w):min(i+w,len));
        
        U = max(thisSubSeq);
        L = min(thisSubSeq);
        
        if (r(i) > U)
            LB = LB + (r(i)-U).^2;
        else
            if (r(i) < L)
                LB = LB + (L-r(i)).^2;
            end
        end
        
    end

end