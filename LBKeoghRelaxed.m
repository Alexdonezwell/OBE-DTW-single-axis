function [LB] = LBKeoghRelaxed(q,r,relativeR,relativeW)

    len = length(q);
    
    w = ceil(len * relativeW);
    rf = ceil(len * relativeR);
    
    LB = 0;   
    for i = rf+1 : len-rf

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