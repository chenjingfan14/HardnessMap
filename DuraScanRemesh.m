%DuraScanRemesh  Refines existing DuraScan mesh with additional points
%
%   Script which returns a ecosWorkflow .spe file and .mat workspace file
%   updated to include new points/mesh identified on the basis of gradient.
%   Saves selected workspace variables to 
%       <Prefix>_<iterations>_setup.mat
%   to the current working directory
%   Saves <Prefix>_<iterations>.spe to a specified subdirectory in the 
%   current working directory, according to <speOut> using writeDuraRows.m
%   Saves a .png file of the last run's results if export_fig is detected ie
%       Output<iterations>.png
%   to the current working directory
%   Requires a *.spe file from the last run
%   Generates 3 figures: the current interpolated results, the current mesh
%   with refinement points targetted and the refined mesh.
%   Requires a workspace with the following variables:
%       D_outline, P_outline - Nx2 arrays of points for the profile outline
%       and datum outline
%       RefIndentDiag,RefIndentHV,RefIndentLoc - corresponding to reference
%       indents made for location purposes
%       distFactor - spacing factor for indents
%       p - current array of points
%       s_pnt - index of p in which the measurements were conducted
%       t - triangular index of p (mesh connectivity)
%       xOff,yOff - x and y value arrays of offset D_outline
%       AddPntsDiag,AddPntsHV,AddPntsHV [optional]: include if there are
%       additional sites on the specimen not contained in RefIndentLoc or p
%   If the specimen has moved since being last run, then a text file with
%   that information is passed via <RefLoc>
%
%   See below for other script variables/parameters.
%
%   Requires xml2struct.m, inpoly.m, gridfit.m, writeDuraRows.m, 
%   rigid_transform_2D.m available on the PATH 
%
%   See also export_fig
%   
%   Copyright 2015 M. J. Roy
%   $Revision: 1.0$  $Date: 2015/10/30$

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Change script variables here
NumRefiningPnts=25;
Prefix='MyPrefix';
speOut=strcat('SomeLocalDirectory\',Prefix); %path & prefix of new spe
LastRunWorkspace='SomeWorkspace.mat';
LastRunResults='SomeLocalDirectory\SomeResults.spe';
RefLoc='MyResultsFile.txt'; %change to valid path/file name accordingly.
%demo=1; %only use for debugging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(LastRunWorkspace);
Iterate=false;
Moved=false;

Results=xml2struct(LastRunResults);

%check if this is the first run since the reference run
%this will be true if there aren't any Result* variables in the loaded
%workspace
if exist('ResultHV','var')==1
    Iterate=true;
end

%check if it's moved
if exist(RefLoc,'file')==2
    Moved=true; l_Datum=dlmread(RefLoc); l_Datum=sortrows(l_Datum,1); 
    l_Datum=l_Datum(:,2:3);
    figure('name','Datum comparison');
    plot(RefIndentLoc(:,1),RefIndentLoc(:,2),'kx'); hold on;
    axis equal
    plot(l_Datum(:,1),l_Datum(:,2),'ro');
    ax=axis;
    %annotation offset on y
    Annot_off_y=(ax(4)-ax(3))*.02;
    for j=1:size(RefIndentLoc,1)
        text(RefIndentLoc(j,1),RefIndentLoc(j,2)-Annot_off_y,3,...
        num2str(j),'FontSize',8,'Color',[0.5 0.5 0.5],...
        'horizontalalignment','center')
    end
    for j=1:size(l_Datum,1)
        text(l_Datum(j,1),l_Datum(j,2)-Annot_off_y,3,...
        num2str(j),'FontSize',8,'Color',[0.5 0.5 0.5],...
        'horizontalalignment','center')
    end
    A=RefIndentLoc;
    B=l_Datum;
    [R,T]=rigid_transform_2D(A,B);
    Aprime=R*A'+repmat(T,1,size(A,1)); Aprime=Aprime';
    %find error in alignment points
    err = sqrt((B(:,1) - Aprime(:,1)).^2+...
        (B(:,2) - Aprime(:,2)).^2);
    [~,ind]=sort(err);
    A=RefIndentLoc(ind(1:3),:);
    B=l_Datum(ind(1:3),:);
    [R,T]=rigid_transform_2D(A,B);
    Aprime=R*A'+repmat(T,1,size(A,1)); Aprime=Aprime';
    err = sqrt((B(:,1) - Aprime(:,1)).^2+...
        (B(:,2) - Aprime(:,2)).^2);
    err = err .* err;
    err = sum(err(:));
    rmse = sqrt(err/size(Aprime,1));
    
    fprintf('RMSE for reorientation: %f mm\n', rmse);
    fprintf('If RMSE is near zero, reorientation was successful.\n');
    plot(Aprime(:,1),Aprime(:,2),'rx');
    legend('Reference location',...
        'New location',...
        'Updated reference location',...
        'location','best')
else
    fprintf('Assuming that the specimen hasn''t moved since the last run\n')
end


% if ~Iterate %then this is the first read after the first run
% get locations, hardness and diagonal
if ~Iterate
    LastRunLoc=zeros(size(p)); %preallocate locations
    LastRunHV=zeros(size(p,1),1); % ... hardness
    LastRunDiag=zeros(size(p,1),1); % ... and diagonals
else
    LastRunLoc=zeros(length(s_pnt{end}),2); %preallocate locations
    LastRunHV=zeros(length(s_pnt{end}),1); % ... hardness
    LastRunDiag=zeros(length(s_pnt{end}),1); % ... and diagonals
end
for j=1:length(Results.Specimen.Row.Point)
    LastRunLoc(j,:)=[str2double(Results.Specimen.Row.Point{j}.XAbs.Text) ...
        str2double(Results.Specimen.Row.Point{j}.YAbs.Text)]./1000;
    LastRunHV(j,1)=str2double...
        (Results.Specimen.Row.Point{j}.Hardness.Text);
    LastRunDiag(j,1)=str2double...
        (Results.Specimen.Row.Point{j}.Diag.Text);
end
Method=Results.Specimen.Row.Point{1}.Method.Text; %for caxis label
if Moved %then translate everything to the new coordinate system
    pPrime=R*p'+repmat(T,1,size(p,1));
    p=pPrime'; 
    LastRunLocPrime=R*LastRunLoc'+repmat(T,1,size(LastRunLoc,1));
    LastRunLoc=LastRunLocPrime'; 
    P_outlinePrime=R*P_outline'+repmat(T,1,size(P_outline,1));
    P_outline=P_outlinePrime';
    D_outlinePrime=R*D_outline'+repmat(T,1,size(D_outline,1));
    D_outline=D_outlinePrime';
    RefIndentLocPrime=R*RefIndentLoc'+repmat(T,1,size(RefIndentLoc,1));
    RefIndentLoc=RefIndentLocPrime'; %for subsequent analyses
    if exist('AddPntsLoc','var')
        AddPntLocPrime=R*AddPntsLoc'+repmat(T,1,size(AddPntsLoc,1));
        AddPntsLoc=AddPntLocPrime'; clear AddPntLocPrime
    end
    clear LastRunLocPrime P_outlinePrime D_outlinePrime RefIndentLocPrime pPrime
end
    
%reorder the Results so that they correspond to the current working mesh.
%overwrite the values for 'p' such that they mirror the actual locations
if ~Iterate
    [~,ind]=sort(s_pnt); %it won't be a cell structure
    p=LastRunLoc(ind,:);
    ResultHV=LastRunHV(ind,:);
    ResultDiag=LastRunDiag(ind,:);
else
    [~,ind]=sort(s_pnt{end}); %it will be cell structure
    p(length(s_pnt{end-1})+1:end,:)=LastRunLoc(ind,:);
    ResultHV=[ResultHV; LastRunHV(ind)];
    ResultDiag=[ResultDiag; LastRunDiag(ind)];
    s_pnt=[s_pnt{1} s_pnt{end}]; %turn the cell into an array
end

%plot the inter/extrapolated result of the last run
p_spread=512;
length_x=max(P_outline(:,1))-min(P_outline(:,1));
length_y=max(P_outline(:,2))-min(P_outline(:,2));
if length_x>length_y
    interp_length=length_x/p_spread;
else
    interp_length=length_y/p_spread;
end
    
x=linspace(min(P_outline(:,1)),max(P_outline(:,1)),round(length_x/interp_length));
y=linspace(min(P_outline(:,2)),max(P_outline(:,2)),round(length_y/interp_length));
if exist('AddPntsLoc','var')==true
    allpnts=[p; RefIndentLoc; AddPntsLoc];
    allHV=[ResultHV;RefIndentHV;AddPntsHV];
else
    allpnts=[p; RefIndentLoc; AddPntsLoc];
    allHV=[ResultHV;RefIndentHV;AddPntsHV];
end
[zg,xg,yg]=gridfit(allpnts(:,1),allpnts(:,2),allHV,x,y);
in_fit=inpoly([reshape(xg,numel(xg),1) reshape(yg,numel(yg),1)],P_outline(:,1:2));
zg(~in_fit)=NaN;

figure('name','Current results');
h=surf(xg,yg,zg); shading flat; hold on;
plot3(p(:,1),p(:,2),ones(size(p(:,1)))+max(ResultHV)+5,'kx');
plot3(P_outline(:,1),P_outline(:,2),ones(size(P_outline(:,1)))+max(ResultHV)+5,'k--',...
    'linewidth',1);
caxis([min(ResultHV) max(ResultHV)]);
view(2); box on; axis equal; h=colorbar('SouthOutside');
% h=colorbar;
xlabel(h,sprintf('%s',Method))
set(gca,'YDir','reverse'); set(gcf,'color','white');


%create a triangulation data structure
TR=triangulation(t,[p ResultHV]); fn=faceNormal(TR);
TR=triangulation(t,p); P=incenter(TR);

%find normals with smallest z components
[~,NormInd]=sortrows(fn,3);
%debug
% plot(P(NormInd(1:NumRefiningPnts),1),P(NormInd(1:NumRefiningPnts),2),...
% 'wo','markeredgecolor','w','markerfacecolor','w'); %debug

%plot all the same stuff as before
%load everything from the _setup file and plot it
theta  = (-pi:pi/12:pi)';

figure('name','Refinement affected elements');
plot(P_outline(:,1),P_outline(:,2),'k--'); hold on; axis equal;
%get axis extents
ax=axis;
%annotation offset on y
Annot_off_y=(ax(4)-ax(3))*.02;
plot(D_outline(:,1),D_outline(:,2),'k-');
for j=1:length(RefIndentDiag)
    circpnts=(RefIndentDiag(j)*distFactor)/2*[cos(theta) sin(theta)];
    patch(circpnts(:,1)+RefIndentLoc(j,1),...
        circpnts(:,2)+RefIndentLoc(j,2),[0.5 0.5 0.5]);
    text(RefIndentLoc(j,1),RefIndentLoc(j,2)-Annot_off_y,3,...
        num2str(j),'FontSize',8,'Color',[0.5 0.5 0.5],...
        'horizontalalignment','center')
end

%plot this
patch('vertices',p,'faces',t,'edgecol','w','facecol',[.8,.9,1]);
%debug
% for j=1:size(RefinementPnts,1)
%     circpnts=(max(RefIndentDiag.*distFactor))/2*[cos(theta) sin(theta)];
%     patch(circpnts(:,1)+RefinementPnts(j,1),...
%         circpnts(:,2)+RefinementPnts(j,2),[0.5 0.5 0.5]);

%     if demo
%         text(p(s_pnt(j),1),p(s_pnt(j),2)-Annot_off_y,3,...
%             num2str(j),'FontSize',8,'Color',[0.1 0.1 0.1],...
%             'horizontalalignment','center')
%     end
% end
%flip y axis
set(gca,'YDir','reverse'); set(gcf,'color','white');

%do the same on the results plot
figure(1)
%debug
% plot3(RefinementPnts(:,1),RefinementPnts(:,2),...
%         ones(size(RefinementPnts(:,1)))+max(ResultHV)+5,'ko',...
%         'markerfacecolor','w','markersize',3);

%{
%Trial of distmesh once refinement has taken place
    bbox=[min(xOff) min(yOff); max(xOff) max(yOff)];
%     fprintf('Now meshing.\nUsing %d seeds on a perimeter of %0.3f mm ...\n', ...
%         nPts, Perimeter);
%     [p,t]=distmesh2d(@dpoly,@huniform,Perimeter/nPts,bbox,...
%         [xOff yOff],[xOff yOff]);
    h0 = mean(sqrt(sum(diff([xOff yOff],1,1).^2,2)));
[newP,newT]=distmesh2d(@dpoly,@huniform,h0*1.5,bbox,[p;RefinementPnts],[xOff yOff]);
%}

figure(2)
%prime while loop
added=true(1,NumRefiningPnts);
refining=true;
i=1;
while refining
    added=true(1,NumRefiningPnts);
        %grab the triangles that need refinement
        if length(NormInd)>=NumRefiningPnts
            RefinementPnts=P(NormInd(1:NumRefiningPnts),:);
        else
            RefinementPnts=P(NormInd,:);
        end
    for j=1:NumRefiningPnts
        %make sure there's enough elements to refine
        try currEl=t(NormInd(j),:);
        catch ME
            fprintf('Ran out of places to put indents.\n')
            refining=false;
            break
        end

        %For each new point, create 3 new elements and one new point to p, required
        %here so that the outline index (O) is respected
        
        %verify that new indents can be placed. Assume that the next indent will be
        %equal to the biggest indent in the element targeted for
        %refinement. If its valid, add it to the end of p, otherwise flag it
        
        addRefinementPnt=false; %like Japanese justice, guilty until innocent
        
        %min hardness in the target element
        local_maxDiag=max(ResultDiag(currEl));
        %create search circle an order of magniture larger than this
        searchPatch=(local_maxDiag*10)/2*[cos(theta) sin(theta)]+...
            repmat(P(NormInd(j),:),length(theta),1);
        %patch(searchPatch(:,1),searchPatch(:,2),[0.85 0.85 0.85]); %debug
        %unique RefinementPnts
        ind=RefinementPnts~=repmat(P(NormInd(j),:),size(RefinementPnts,1),1);
        if exist('AddPntsLoc','var')==true
            testedPoints=[p;...
                RefinementPnts(ind(:,1)&ind(:,2),:);...
                RefIndentLoc;...
                AddPntsLoc];
        else
            testedPoints=[p;...
                RefinementPnts(ind(:,1)&ind(:,2),:);...
                RefIndentLoc];
        end
        testedPoints=[p;RefinementPnts(ind(:,1)&ind(:,2),:);RefIndentLoc];
        in=inpoly(testedPoints,...
            searchPatch); %all potential points
        local_p=testedPoints(in,:);
        if ~isempty(local_p)
            currDists=zeros(size(local_p,1)); %preallocate
            for k=size(local_p,1)
                %find the linear distance to all of the points in the search
                %patch
                currDists(k)=ldist([P(NormInd(j),:) local_p(k,:)]);
            end
            if currDists>distFactor*local_maxDiag*1 %increase from 1 to make conservative
                addRefinementPnt=true;
            end
        else addRefinementPnt=true;
        end
        if addRefinementPnt
            p(end+1,:)=P(NormInd(j),:);
            ThesePointsGetAdded(i,:)=P(NormInd(j),:);
            i=i+1;
            t=[t;
                currEl(1) currEl(2) size(p,1);
                currEl(2) currEl(3) size(p,1);
                currEl(3) currEl(1) size(p,1)];
        else
            added(j)=0;
        end
    end
    if refining
        fprintf('%d out of %d refinement points added ...\n',nnz(added),NumRefiningPnts);
        if nnz(added)~=NumRefiningPnts
            %then remove the tested entries from NormInd
            NormInd(1:NumRefiningPnts)=[];
            %refine NumRefiningPnts
            NumRefiningPnts=NumRefiningPnts-nnz(added);
        else
            break
        end
    end

end %while true

% Find all edges in mesh, note internal edges are repeated
F=TR.ConnectivityList;
E = sort([F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)]')';
% determine uniqueness of edges
[u,m,n] = unique(E,'rows');
% determine counts for each unique edge
counts = accumarray(n(:), 1);
% extract edges that only occurred once
O = u(counts==1,:);
dt = delaunayTriangulation(p(:,1),p(:,2),O);

%now minimize distance between the refinement points
%create local entries to s_pnt
os_pnt=length(s_pnt); %original length of s_pnt
new_s_pnt=NaN(1,size(ThesePointsGetAdded,1));

dist=pdist2(ThesePointsGetAdded,ThesePointsGetAdded,'cityblock');
dist2=dist;
new_s_pnt(1)=1; %first point will be the first RefinementPnt
for j=2:size(ThesePointsGetAdded,1)
    dist(:,new_s_pnt(j-1))=Inf;
    [~,closest_idx]=min(dist(new_s_pnt(j-1),:));
    new_s_pnt(j)=closest_idx;
    %debug as follows
%     plot([p(new_s_pnt(j-1)+os_pnt,1) p(new_s_pnt(j)+os_pnt,1) ],...
%         [p(new_s_pnt(j-1)+os_pnt,2) p(new_s_pnt(j)+os_pnt,2)],'r-');
%     pause(0.1);
%     text(p(new_s_pnt(j)+os_pnt,1),p(new_s_pnt(j)+os_pnt,2)-Annot_off_y,3,...
%             num2str(j),'FontSize',8,'Color',[0.1 0.1 0.1],...
%             'horizontalalignment','center')
end


for j=1:size(p,1)
    circpnts=(max(RefIndentDiag)*distFactor)/2*[cos(theta) sin(theta)];
    patch(circpnts(:,1)+p(j,1),...
        circpnts(:,2)+p(j,2),[0.5 0.5 0.5]);
    %debug
%     if demo
%         text(p(s_pnt(j),1),p(s_pnt(j),2)-Annot_off_y,3,...
%             num2str(j),'FontSize',8,'Color',[0.1 0.1 0.1],...
%             'horizontalalignment','center')
%     end
end


figure('name','Refined mesh');
plot(P_outline(:,1),P_outline(:,2),'k--'); hold on;
IO=isInterior(dt); t=dt.ConnectivityList(IO,:);
plottedMesh=patch('vertices',p,'faces',t,'edgecol','w','facecol',[.8,.9,1]);
plot([p(O(:,1),1) p(O(:,2),1)]',[p(O(:,1),2) p(O(:,2),2)]','k-');
for j=1:length(s_pnt)
    circpnts=(ResultDiag(j).*distFactor)/2*[cos(theta) sin(theta)];
    patch(circpnts(:,1)+p(j,1),...
        circpnts(:,2)+p(j,2),[0.5 0.5 0.5]);
end

for j=1:length(new_s_pnt)
    circpnts=(min(RefIndentDiag)*distFactor)/2*[cos(theta) sin(theta)];
    patch(circpnts(:,1)+p(new_s_pnt(j)+os_pnt,1),...
        circpnts(:,2)+p(new_s_pnt(j)+os_pnt,2),[1 0.5 0.5]);
end
plot(RefinementPnts(~added,1),RefinementPnts(~added,2),'kx')

axis equal;
set(gca,'YDir','reverse'); set(gcf,'color','white');


%now write a new programme with these new points
if ~Iterate
    iterations=1;
else
    iterations=iterations+1;
end
writeDuraRows(p(new_s_pnt+os_pnt,:),...
    strcat(speOut,'_',num2str(iterations),'.spe'));

%write a new workspace, dealing with multiplicity via cells
s_pnt={s_pnt new_s_pnt};

%handle additional points
if exist('AddPntsLoc','var')==true
    save(strcat(Prefix,'_',num2str(iterations),'_setup.mat'),...
        'RefIndentLoc', 'RefIndentHV','RefIndentDiag',...
        'P_outline','D_outline','xOff','yOff','p','t','s_pnt',...
        'ResultHV','ResultDiag','distFactor','iterations',...
        'AddPntsDiag','AddPntsHV','AddPntsLoc');
else
    save(strcat(Prefix,'_',num2str(iterations),'_setup.mat'),...
        'RefIndentLoc', 'RefIndentHV','RefIndentDiag',...
        'P_outline','D_outline','xOff','yOff','p','t','s_pnt',...
        'ResultHV','ResultDiag','distFactor','iterations');
end


if exist('export_fig','file')==2
    if Moved
        figure(2);
    else
        figure(1);
    end
    export_fig(sprintf('Output%d.png',iterations-1),'-a2','-m2')
end
