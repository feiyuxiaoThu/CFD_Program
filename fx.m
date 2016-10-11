function [ ut0 ] = fx( x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if x<0.3
    ut0 = 0;
elseif  (x<=0.6) && (x>=0.3)
     ut0 = 1;
else
        ut0 = 1 + 2.5* (x - 0.6);
end

        
end

    


