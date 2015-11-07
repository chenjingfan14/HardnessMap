%respace_equally  respace points equally around a given 2D hull
%   [X,Y,P,Q]=respace_equally(B,I) will create Qx1 arrays X and Y
%   of the Nx2 input point matrix B of evenly spaced points. The
%   number of points Q is returned along with the perimeter, P.
%   Input I is either an int8 variable specifying the number of points
%   ie Q=I, or a target length between [X1 YI],[X2 Y2] etc.
%
%   Copyright 2015 M. J. Roy
%   $Revision: 1.0$  $Date: 2015/10/30$
function [xNew,yNew,Perimeter,nPts]=respace_equally(B,input)
%Get total distance between each point and then take cumulative sum;
%parametric coordinate that makes small steps for points close together and
%larger for ones further apart.
distance = sqrt(sum(diff(B,1,1).^2,2));  %# Distance between subsequent points
s = [0; cumsum(distance)];               %# Parametric coordinate
%Solve for the perimeter
Perimeter=sum(distance);
if ~isinteger(input)
    nPts=round(Perimeter/input);
else
    nPts=input;
end

%interpolate for a new set of points equally spaces along lines joining
%points

sNew = linspace(0,s(end),nPts).';   %'# nPts evenly spaced points from 0 to s(end)
xNew = interp1q(s,B(:,1),sNew);     %# Interpolate new x values
yNew = interp1q(s,B(:,2),sNew);     %# Interpolate new y values
nPts=double(nPts);
