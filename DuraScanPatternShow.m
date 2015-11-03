%DuraScanPatternShow  Interpolates hardness results from a hardness map
%
%   Script which returns an image which can then be output via export_fig
%   Saves selected workspace variables to intermediate_<Prefix>.mat to the
%   current working directory
%   Saves a .png file of the results if export_fig is detected.
%   Requires a <result> .spe file from the DuraScan
%   Requires an <outline> file which contains Nx2 points, whitespace
%   delimited, one point per line.
%   See below for other script variables/parameters.
%
%   Requires xml2struct.m, inpoly.m, gridfit.m available on the PATH 
%
%   See also export_fig.
%   
%   Copyright 2015 M. J. Roy
%   $Revision: 1.0$  $Date: 2015/10/30$

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Change script variables here
R=0; %rotation of the outline/results about the centroid (degrees)
outline='SomeDirectory\MyOutline.txt'; %path and file name for outline
result='SomeDirectory\MyResults.spe'; %path and file name for corresponding *.spe
outprefix='SomePrefix'; %Prefix of output intermediate_<Prefix>.mat
redo = 1; %boolean - true: discard existing .mat, -false: load .mat file
ScaleBarLength=5; %mm (or 0 for no outline)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%check whether to load .mat file or not; 
%delete .mat file if redo is true and the file exists
if redo && exist(strcat('intermediate_',outprefix,'.mat'),'file') ==2 
    delete(strcat('intermediate_',outprefix,'.mat'));
end

%check if the .mat file exists, if it doesn't, start processing
if ~exist(strcat('intermediate_',outprefix,'.mat'),'file')
    A=xml2struct(result); %call xml2struct
    Points=A.Specimen.Row.Point; %store all 'point' data in a cell structure

    %For all entries in Points, extract the hardness and locations
    P_abs=zeros(length(Points),2); Hardness=zeros(length(Points),1);
    for j=1:length(Points)
        P_abs(j,:)=[str2double(Points{j}.XAbs.Text) ...
            str2double(Points{j}.YAbs.Text)]./1000;
        Hardness(j,1)=str2double(A.Specimen.Row.Point{j}.Hardness.Text);
    end
    Method=Points{1}.Method.Text;
else
    load(strcat('intermediate_',outprefix,'.mat'))
end %fi intermediate exists

O_pnt=dlmread(outline); %read in the outline
O_pnt(end+1,:)=O_pnt(1,:); %close the outline

Rm=[cosd(R) -sind(R); sind(R) cosd(R)]; %create rotation matrix

cent=mean(O_pnt); %calculate centroid

%translate all points to the centroid of the outline
O_pnt=[O_pnt(:,1)-cent(1) O_pnt(:,2)-cent(2)]; 
P_abs=[P_abs(:,1)-cent(1) P_abs(:,2)-cent(2)];

%rotate them about the z axis by R
O_pnt=O_pnt*Rm;
P_abs=P_abs*Rm;


%set up interpolation vectors
p_spread=512; %max number of points in longest direction
length_x=max(O_pnt(:,1))-min(O_pnt(:,1));
length_y=max(O_pnt(:,2))-min(O_pnt(:,2));
if length_x>length_y
    interp_length=length_x/p_spread;
else
    interp_length=length_y/p_spread;
end

%interpolation vectors
x=linspace(min(O_pnt(:,1)),max(O_pnt(:,1)),round(length_x/interp_length));
y=linspace(min(O_pnt(:,2)),max(O_pnt(:,2)),round(length_y/interp_length));

%acquire interpolated values from gridfit
[zg,xg,yg]=gridfit(P_abs(:,1),P_abs(:,2),Hardness,x,y);
%find values which are within the outline
in_fit=inpoly(...
    [reshape(xg,numel(xg),1) reshape(yg,numel(yg),1)],O_pnt(:,1:2));
%kill points outside for surface plotting
zg(~in_fit)=NaN;

%plot stuff on fullscreen figure
figure('units','normalized','outerposition',[0 0 1 1],...
    'name',sprintf('%s results',outprefix));
surf(xg,yg,zg); shading flat; hold on;
%bump up annotations on z axis to prevent them being overwritten
plot3(P_abs(:,1),P_abs(:,2),ones(size(P_abs(:,1)))+max(Hardness)+5,'kx');
plot3(O_pnt(:,1),O_pnt(:,2),ones(size(O_pnt(:,1)))+max(Hardness)+5,'k-',...
    'linewidth',1.5);
%set the colorbar limits
caxis([min(Hardness) max(Hardness)]);
%turn off the axis, flip the y direction and make the background white
set(gca,'YDir','reverse'); %matches orientation on the tester
axis off;
set(gcf,'color','white');

%Create the scalebar
if ScaleBarLength~=0
    ax=axis;
    Lo=0.05.*[ax(2)-ax(1) ax(4)-ax(3)];
    To=[-Lo(1) -Lo(2);
    Lo(1) -Lo(2);
    -Lo(1) Lo(2);
    Lo(1) Lo(2)];

    ScaleB=plot([max(O_pnt(:,1))-Lo(1)...
        max(O_pnt(:,1))-Lo(1)-ScaleBarLength],...
          [min(O_pnt(:,2))-Lo(2) min(O_pnt(:,2))-Lo(2)],'-','LineWidth',3);
      set(ScaleB,'color','r');
      text(max(O_pnt(:,1))-Lo(1)-0.5*ScaleBarLength,...
          min(O_pnt(:,2))-1.5*Lo(2),...
          sprintf('%d mm',ScaleBarLength),...
          'color','k','fontsize',12,...
          'horizontalalignment','center','Interpreter','none');

end
%put a label on the results
h=colorbar;
xlabel(h,sprintf('%s',Method))
axis equal;
view(2);

if exist('export_fig','file')==2
    export_fig(sprintf('%s',strcat(outprefix,'.png')),'-a2','-m2');
end