function [r] = ranGenerator( lb, up, n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 3     
        r = (up - lb).*rand([1,1]) + lb;
    else
        r = (up - lb).*rand([1,n]) + lb;
    end
end

