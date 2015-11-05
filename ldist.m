%ldist  calculates the length of a 2D line segment 
%   L=ldist(LINE) where LINE=[x1 y1 x2 y2]
%
%   Copyright 2015 M. J. Roy
%   $Revision: 1.0$  $Date: 2015/10/30$

function L=ldist(l)
dx=l(3)-l(1);
dy=l(4)-l(2);
L=sqrt(dx^2+dy^2);
