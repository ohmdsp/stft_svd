function [b,element]=splitoff(a,position)

% This function separates the top or bottom element of input vector and places
% remainder into output vector.
%
%     [b,element]=splitoff(a,position)                               version 1.0
%
% a        -- input vector
% position -- set to 'top' to split off top element, or 'bottom' to split off
%                  bottom element
% b        -- truncated version of input vector
% element  -- top or bottom element of input vector

n = length(a);
if strcmp(position,'top')
    element = a(1);
    b = a(2:n);
elseif strcmp(position,'bottom')
    element = a(n);
    b = a(1:n-1);
end
%
