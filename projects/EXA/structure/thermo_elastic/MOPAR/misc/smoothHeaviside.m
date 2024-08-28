function [y,dy,dy2]=smoothHeaviside(x,thres,beta)
%
% smoothed Heaviside funtion
%
% Input: x     - independent variable
%        thres - threshold
%        beta  - strengh of projection
%

if beta < 0
    y=x;   dy=1; dy2=0; return;
end

y   = (tanh(beta*thres)+tanh(beta*(x - thres)))/(tanh(beta*thres)+tanh(beta*(1 - thres)));
dy  = (beta*(tanh(beta*(thres - x))^2 - 1))/(tanh(beta*(thres - 1)) - tanh(beta*thres));
dy2 = (2*beta^2*tanh(beta*(thres - x))*(tanh(beta*(thres - x))^2 - 1))/(tanh(beta*(thres - 1)) - tanh(beta*thres));
