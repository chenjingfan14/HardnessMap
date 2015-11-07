function q=simpqual(p,t,type)
%SIMPQUAL Simplex quality.
%   Q=SIMPQUAL(P,T,TYPE)
%
%   TYPE == 1: Radius Ratio (default)
%   TYPE == 2: Approximate

%   Copyright (C) 2004-2012 Per-Olof Persson.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT
% DistMesh is a collection of MATLAB functions for generation and
% manipulation of unstructured meshes. DistMesh is Copyright (C) 2004-2012
% Per-Olof Persson.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

% If you use DistMesh in any program or publication, please acknowledge
% its authors by adding a reference to: Per-Olof Persson and Gilbert
% Strang, "A Simple Mesh Generator in MATLAB," SIAM Review Vol. 46 (2)
% 2004.
%-------------------------------------------------------- 
% Modified by M. J. Roy, University of Manchester, 15-5-2015
% For use in DuraScan programming suite
% Fair use under GNU license.

if nargin<3
  type=1;
end

switch type
 case 1
  % RADIUS RATIO
  switch size(p,2)
   case 1
    q=ones(1,size(t,2));
   case 2
    a=sqrt(sum((p(t(:,2),:)-p(t(:,1),:)).^2,2));
    b=sqrt(sum((p(t(:,3),:)-p(t(:,1),:)).^2,2));
    c=sqrt(sum((p(t(:,3),:)-p(t(:,2),:)).^2,2));
    r=1/2*sqrt((b+c-a).*(c+a-b).*(a+b-c)./(a+b+c));
    R=a.*b.*c./sqrt((a+b+c).*(b+c-a).*(c+a-b).*(a+b-c));
    q=2*r./R;
   case 3
    d12=p(t(:,2),:)-p(t(:,1),:);
    d13=p(t(:,3),:)-p(t(:,1),:);
    d14=p(t(:,4),:)-p(t(:,1),:);
    d23=p(t(:,3),:)-p(t(:,2),:);
    d24=p(t(:,4),:)-p(t(:,2),:);
    d34=p(t(:,4),:)-p(t(:,3),:);
    v=abs(dot(cross(d12,d13,2),d14,2))/6;
    s1=sqrt(sum(cross(d12,d13,2).^2,2))/2;
    s2=sqrt(sum(cross(d12,d14,2).^2,2))/2;
    s3=sqrt(sum(cross(d13,d14,2).^2,2))/2;
    s4=sqrt(sum(cross(d23,d24,2).^2,2))/2;
    p1=sqrt(sum(d12.^2,2)).*sqrt(sum(d34.^2,2));
    p2=sqrt(sum(d23.^2,2)).*sqrt(sum(d14.^2,2));
    p3=sqrt(sum(d13.^2,2)).*sqrt(sum(d24.^2,2));
    q=216*v.^2./(s1+s2+s3+s4)./sqrt((p1+p2+p3).*(p1+p2-p3).* ...
                                    (p1+p3-p2).*(p2+p3-p1));
   otherwise
    error('Dimension not implemented.');
  end
 case 2
  % APPROXIMATE FORMULA
  switch size(p,2)
   case 1
    q=ones(1,size(t,2));
   case 2
    d12=sum((p(t(:,2),:)-p(t(:,1),:)).^2,2);
    d13=sum((p(t(:,3),:)-p(t(:,1),:)).^2,2);
    d23=sum((p(t(:,3),:)-p(t(:,2),:)).^2,2);
    q=4*sqrt(3)*abs(simpvol(p,t))./(d12+d13+d23);
   case 3
    d12=sum((p(t(:,2),:)-p(t(:,1),:)).^2,2);
    d13=sum((p(t(:,3),:)-p(t(:,1),:)).^2,2);
    d14=sum((p(t(:,4),:)-p(t(:,1),:)).^2,2);
    d23=sum((p(t(:,3),:)-p(t(:,2),:)).^2,2);
    d24=sum((p(t(:,4),:)-p(t(:,2),:)).^2,2);
    d34=sum((p(t(:,4),:)-p(t(:,3),:)).^2,2);
    q=216*abs(simpvol(p,t))/sqrt(3)./(d12+d13+d14+d23+d24+d34).^(3/2);
   otherwise
    error('Dimension not implemented.');
  end
 otherwise
  error('Incorrect type.');
end
