function [ u1t  ] = bt( chioce,t )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    if chioce == 1
        u1t = 2;
    end
    
    if chioce == 2
        u1t = 2 + sin(10*t);
            
    end

end

