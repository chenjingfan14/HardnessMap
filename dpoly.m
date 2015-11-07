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
function d=dpoly(p,pv)


np=size(p,1);
nvs=size(pv,1)-1;

ds=dsegment(p,pv);
%ds=zeros(np,nvs);
%for iv=1:nvs
%  ds(:,iv)=donesegment(p,pv(iv:iv+1,:));
%end
d=min(ds,[],2);

d=(-1).^(inpolygon(p(:,1),p(:,2),pv(:,1),pv(:,2))).*d;

% MEXED

%function ds=donesegment(p,pv)
%
%e=ones(size(p,1),1);
%
%v=diff(pv,1);
%w=p-e*pv(1,:);
%
%c1=sum(w.*v(e,:),2);
%c2=sum(v(e,:).^2,2);
%
%ds=0*e;
%
%ix=c1<=0;
%ds(ix)=sqrt(sum((p(ix,:)-pv(1*ones(sum(ix),1),:)).^2,2));
%
%ix=c1>=c2;
%ds(ix)=sqrt(sum((p(ix,:)-pv(2*ones(sum(ix),1),:)).^2,2));
%
%ix=c1>0 & c2>c1;
%nix=sum(ix);
%if nix>0
%  Pb=ones(nix,1)*pv(1,:)+c1(ix)./c2(ix)*v;
%  ds(ix)=sqrt(sum((p(ix,:)-Pb).^2,2));
%end
