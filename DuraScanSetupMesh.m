%DuraScanSetupMesh  Initializes a hardness map
%
%   Script which returns a ecosWorkflow .spe file and .mat workspace file
%   which includes experiment-specific points and a mesh.
%   Saves selected workspace variables to 
%       <Prefix>_setup.mat
%   to the current working directory
%   Saves <Prefix>.spe to a specified subdirectory in the 
%   current working directory, according to <speOut> using writeDuraRows.m
%   Generates 2 figures: the mesh as it is calculated, and the final meshed
%   instance with all current indents shown.
%   Generates a workspace with the following variables:
%       D_outline, P_outline - Nx2 arrays of points for the profile outline
%       and domain outline
%       RefIndentDiag,RefIndentHV,RefIndentLoc - corresponding to reference
%       indents made for location purposes
%       distFactor - spacing factor for indents
%       p - current array of points
%       s_pnt - index of p in which the measurements will be conducted
%       t - triangular index of p (mesh connectivity)
%       xOff,yOff - x and y value arrays of offset D_outline
%       AddPntsDiag,AddPntsHV,AddPntsHV [if read]: include if there are
%       additional sites on the specimen not contained in RefIndentLoc or p
%   Requires a reference .spe file containing the results of the setup
%   indents
%   Requires a perimeter/profile file in .txt format, whitespace delimited,
%   one point per row. Order is not important.
%   Requires a local (or domain) outline in the same format as the
%   perimeter. This outline will be meshed.
%   
%   See below for other script variables/parameters.
%
%   Requires writeDuraRows.m, clipper.mexw64, xml2struct.m, 
%   respace_equally.m, inpoly.m, distmesh2D (and its dependencies)
%       -dpoly.m
%       -huniform.m
%       -simpqual.m
%       -dsegment.mexw64
%
%   See also distmesh2D
%   
%   Copyright 2015 M. J. Roy
%   $Revision: 1.0$  $Date: 2015/10/30$

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Change script variables here
distFactor=3; %default
seed=2; %mm, or set seed=int8(number of seeds)
Prefix='DemoReprocessor';
speOut=strcat('Programs\',Prefix); %path & prefix of new spe
RefIndent='Results_Outlines\TryReorient.spe';
PeriOutline='Results_Outlines\OverviewOutline_Conv.txt';
DomainOutline='Results_Outlines\EvaluationOutline_Conv.txt';
AddPoints='Results_Outlines\Already_done.txt';
AddPointsLoad=1;
demo=0; %make equal to 1 or true to see indent sequence
doOver=1; %if true, it will only load the contents of <Prefix>_setup.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if doOver && exist(strcat(Prefix,'_setup.mat'),'file')==2
    delete(strcat(Prefix,'_setup.mat'));
end

if exist(strcat(Prefix,'_setup.mat'),'file')~=2

    %read in profiles
    P_outline=dlmread(PeriOutline);
    D_outline=dlmread(DomainOutline);
    
    %check for additional points
    ap=false;%flag for subsequent analysis
    if exist(AddPoints,'file')==2 %if a valid file is found
        fprintf('Accounting for existing indents beyond reference ...\n');
        ap=true; 
        %read them
        AddPnts=dlmread(AddPoints);

        AddPntsDiag=sqrt(1.8544*AddPointsLoad./(AddPoints(:,end)));
        AddPntsLoc=AddPoints(:,1:2);
        AddPntsHV=AddPoints(:,end);
    end
    
    RawRefIndent=xml2struct(RefIndent);
    %get locations, hardness and diagonal
    for j=1:length(RawRefIndent.Specimen.Row.Point)
        RefIndentLoc(j,:)=[str2double(RawRefIndent.Specimen.Row.Point{j}.XAbs.Text) ...
            str2double(RawRefIndent.Specimen.Row.Point{j}.YAbs.Text)]./1000;
        RefIndentHV(j,1)=str2double...
            (RawRefIndent.Specimen.Row.Point{j}.Hardness.Text);
        RefIndentDiag(j,1)=str2double...
            (RawRefIndent.Specimen.Row.Point{j}.Diag.Text);
    end
    %for subsequent plotting
    theta  = (-pi:pi/12:pi)';
    
    figure; 
    set(gca,'YDir','reverse'); %to match durascan orientation
    plot(P_outline(:,1),P_outline(:,2),'k--'); hold on; axis equal;
    %get axis extents
    ax=axis;
    %annotation offset on y
    Annot_off_y=(ax(4)-ax(3))*.02;
    %plot the domain outline
    plot(D_outline(:,1),D_outline(:,2),'k-');
    
    %plot the reference points and their order
    for j=1:length(RefIndentDiag)
        %plastically affected zone is circle with radius (diag*distfactor)/2
        circpnts=(RefIndentDiag(j)*distFactor)/2*[cos(theta) sin(theta)];
        patch(circpnts(:,1)+RefIndentLoc(j,1),...
            circpnts(:,2)+RefIndentLoc(j,2),[0.5 0.5 0.5]);
        text(RefIndentLoc(j,1),RefIndentLoc(j,2)-Annot_off_y,3,...
            num2str(j),'FontSize',8,'Color',[0.5 0.5 0.5],...
            'horizontalalignment','center')
    end
    
    %plot any additional points
    if ap
        for j=1:length(AddPntsDiag)
            %plastically affected zone is circle with radius (diag*distfactor)/2
            circpnts=(AddPntsDiag(j)*distFactor)/2*[cos(theta) sin(theta)];
            patch(circpnts(:,1)+AddPntsLoc(j,1),...
                circpnts(:,2)+AddPntsLoc(j,2),[0.5 0.5 0.5]);
            text(AddPntsLoc(j,1),AddPntsLoc(j,2)+Annot_off_y,3,...
                num2str(j),'FontSize',8,'Color',[1 0.5 0.5],...
                'horizontalalignment','center');
        end
    end %fi ap
    
    %have to convert to int64 for clipper, use an arbitrary scaling factor
    scale=2^15;
    %offset the outline by the largest indent diagonal in the reference
    %points
    if ~ap %no additional points
        Offset=max(RefIndentDiag.*distFactor);
    else
        Offset=max([RefIndentDiag;AddPntsDiag].*distFactor);
    end
    Inpol.x=int64(D_outline(:,1)*scale); Inpol.y=int64(D_outline(:,2)*scale);
    Outpol=clipper(Inpol,-Offset*scale,1);

    %convert Outpol back to real
    xOff=[Outpol(1).x; Outpol(1).x(1)]/scale;
    yOff=[Outpol(1).y; Outpol(1).y(1)]/scale;
    
    %respace the domain for meshing
    [xOff,yOff,Perimeter,nPts]=respace_equally([xOff yOff],seed);
    
    %set up mesh stats
    fstats=@(p,t) fprintf('%d nodes, %d elements, min quality %.2f\n', ...
        size(p,1),size(t,1),min(simpqual(p,t)));
    %get bounding box
    bbox=[min(xOff) min(yOff); max(xOff) max(yOff)];
    %mesh it
    fprintf('Now meshing.\nUsing %d seeds on a perimeter of %0.3f mm ...\n', ...
        nPts, Perimeter);
    [p,t]=distmesh2d(@dpoly,@huniform,Perimeter/nPts,bbox,...
        [xOff yOff],[xOff yOff]);
    
    fstats(p,t);
    figure(1)

    %check if reference/additional points are within where further 
    %points are intended
    if ap
        searchPnts=[AddPntsLoc; RefIndentLoc];
        searchDiag=[AddPntsDiag; RefIndentDiag];
    else
        searchPnts=RefIndentLoc; searchDiag=RefIndentDiag;
    end
    %calculate the distance and index between all points
    [k,RefDist] = dsearchn(p,t,searchPnts);
    DeletePointTrue=true(size(k)); %initialization
    for j=1:length(k)
        %if there are any points too close
        if RefDist(j)<=searchDiag*distFactor 
            t(any(t==k(j),2),:)=[]; %delete their index in the mesh
        else
            DeletePointTrue(j)=0; %adjust the deletepoint entry
        end
    end
    
    %for all points in the mesh that are too close to existing points,
    %delete them
    p(k(DeletePointTrue),:)=[]; RefDist(DeletePointTrue)=[];
    
    %make sure that the number of values in p is not exceeded by t's index
    t=delaunayn(p); %calculate delaunay
    pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
    in=inpoly(pmid,[xOff yOff]);
    t=t(in,:);
    %plot this updated mesh
    plottedMesh=patch('vertices',p,'faces',t,'edgecol','w','facecol',[.8,.9,1]);
    alpha(plottedMesh,0.95);
    
    N=size(p,1);
    %Minimize total path through points with pdist2
    dist=pdist2(p,p,'cityblock'); dist2=dist;
    %check if there's any positions that invalidate hardness
    dist2(dist2==0)=NaN;
    minDist=min(dist2,[],1);
    fprintf...
        ('Min. dist. between points: %0.3f\n',min([minDist RefDist']));
    fprintf('Min. based on reference values: %0.3f\n',max(RefIndentDiag.*distFactor));
    fprintf('Largest distance: %0.3f\n',max(minDist));
    fprintf('St. dev. of dist: %0.3f\n',std(minDist));
    
    
    s_pnt=NaN(1,N);
    s_pnt(1)=1; %start with the first entry in p
    
    %pathdist=0; % for total path length between points
    for j=2:N
        dist(:,s_pnt(j-1))=Inf;
        [~,closest_idx]=min(dist(s_pnt(j-1),:));
        s_pnt(j)=closest_idx;
    end
    
   for j=1:size(p,1)
       %plastically affected zone is circle with radius diag*distfactor/2
       circpnts=(max(searchDiag)*distFactor)/2*[cos(theta) sin(theta)];
       %draw them
       patch(circpnts(:,1)+p(j,1),...
        circpnts(:,2)+p(j,2),[0.5 0.5 0.5]);
    if demo %display the order of measurements?
        text(p(s_pnt(j),1),p(s_pnt(j),2)-Annot_off_y,3,...
            num2str(j),'FontSize',8,'Color',[0.1 0.1 0.1],...
            'horizontalalignment','center')
    end
   end
   set(gca,'YDir','reverse'); set(gcf,'color','white');
    %hardness measurement points are now according to s_pnt which means that
    %H(1) corresponds to p(s_pnt(1)), H(2) will be p(s_pnt(2)) etc.

    %save these points
    if ap
        save(strcat(Prefix,'_setup.mat'),...
            'RefIndentLoc', 'RefIndentHV','RefIndentDiag',...
            'AddPntsLoc','AddPntsHV','AddPntsDiag',...
            'P_outline','D_outline','xOff','yOff','p','t',...
            's_pnt','distFactor');
    else
        save(strcat(Prefix,'_setup.mat'),...
            'RefIndentLoc', 'RefIndentHV','RefIndentDiag',...
            'P_outline','D_outline','xOff','yOff','p','t',...
            's_pnt','distFactor');
    end
else
    %load everything from the _setup file and plot it
    theta  = (-pi:pi/12:pi)';
    load(strcat(Prefix,'_setup.mat'));
    figure; 
    plot(P_outline(:,1),P_outline(:,2),'k--'); hold on; axis equal;
    %get axis extents
    ax=axis;
    %annotation offset on y
    Annot_off_y=(ax(4)-ax(3))*.02;
    plot(D_outline(:,1),D_outline(:,2),'k-');
    %plot reference locations
    for j=1:length(RefIndentDiag)
        %plastically affected zone is circle with radius (diag*distfactor)/2
        circpnts=(RefIndentDiag(j)*distFactor)/2*[cos(theta) sin(theta)];
        patch(circpnts(:,1)+RefIndentLoc(j,1),...
            circpnts(:,2)+RefIndentLoc(j,2),[0.5 0.5 0.5]);
        text(RefIndentLoc(j,1),RefIndentLoc(j,2)-Annot_off_y,3,...
            num2str(j),'FontSize',8,'Color',[0.5 0.5 0.5],...
            'horizontalalignment','center')
    end
    
    %plot any additional points
    if ap
        for j=1:length(AddPntsDiag)
            %plastically affected zone is circle with radius (diag*distfactor)/2
            circpnts=(AddPntsDiag(j)*distFactor)/2*[cos(theta) sin(theta)];
            patch(circpnts(:,1)+AddPntsLoc(j,1),...
                circpnts(:,2)+AddPntsLoc(j,2),[0.5 0.5 0.5]);
            text(AddPntsLoc(j,1),AddPntsLoc(j,2)+Annot_off_y,3,...
                num2str(j),'FontSize',8,'Color',[1 0.5 0.5],...
                'horizontalalignment','center');
        end
    end %fi ap

    plottedMesh=patch('vertices',p,'faces',t,'edgecol','w','facecol',[.8,.9,1]);
    for j=1:size(p,1)
        %plastically affected zone is circle with radius (diag*distfactor)/2
        circpnts=(max(RefIndentDiag)*distFactor)/2*[cos(theta) sin(theta)];
        patch(circpnts(:,1)+p(j,1),...
            circpnts(:,2)+p(j,2),[0.5 0.5 0.5]);
        if demo
            text(p(s_pnt(j),1),p(s_pnt(j),2)-Annot_off_y,3,...
                num2str(j),'FontSize',8,'Color',[0.1 0.1 0.1],...
                'horizontalalignment','center')
        end
    end
    set(gca,'YDir','reverse'); set(gcf,'color','white');
end %if the setup file exists and redo isn't true

%write a programme
writeDuraRows(p(s_pnt,:),...
    strcat(speOut,'.spe'));
