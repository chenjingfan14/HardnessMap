%writeDuraRows  writes an .spe file to be read into ecosWorkFlow
%   writeDuraRows(P,F) will write the Nx2 matrix of x & y coordinates
%   to file targeted with string F with defaults HV 1/20x objective.
%
%   writeDuraRows(P,F,M) will do the same as the above with string M
%   representing the method (ie 'HV 1', 'HV 0,5' etc).
%
%   writeDuraRows(P,F,M,O) will do the same as the above with string O
%   representing the objective (ie '20x', '40x' etc).
%
%   
%   Copyright 2015 M. J. Roy
%   $Revision: 1.0$  $Date: 2015/10/30$

function writeDuraRows(varargin)
if length(varargin)<2 || length(varargin)>5
    fprintf('Inputs for writeDuraRows not specified correctly.')
    return;
end
p=varargin{1};
fname=varargin{2};
if length(varargin)>2
    HVstr=varargin{3};
else
    HVstr='HV 1';
end

if length(varargin)>3
    Objstr=varargin{4};
else
    Objstr='20x';
end

Ipath=fullfile(pwd,strcat('Title',datestr(now),'.jpg'));
User='User';


AbsMeasPnts=round(p.*1000); %input as mm, output as microns
GlobalStartPoint=round(mean(p).*1000);
RelMeasPnts=...
    [AbsMeasPnts(:,1)-GlobalStartPoint(1) AbsMeasPnts(:,2)-GlobalStartPoint(2)];
RelMeasPnts=RelMeasPnts./1000;

dt=datestr(now,'mm/dd/yyyy HH:MM:SS AM');


zoom=1;
Commentstr='Enter Comment (Optional)';


%Uses fprintf to write an .spe file (.xml) instead of xml protocols Matlab
%ships with are insufficient
fid=fopen(fname,'w+');
fprintf(fid,'<?xml version="1.0"?>\r\n'); %header
fprintf(fid,'<Specimen>\r\n'); %start specimen field

fprintf(fid,'  <Testtype>Series Measurement</Testtype>\r\n'); %type of test

fprintf(fid,'  <OCImagePath>%s</OCImagePath>\r\n',Ipath); %set path

fprintf(fid,'  <SpecimenStartPoint>\r\n'); %do global start point
fprintf(fid,'    <XAbs>%d</XAbs>\r\n',GlobalStartPoint(1));
fprintf(fid,'    <YAbs>%d</YAbs>\r\n',GlobalStartPoint(2));
fprintf(fid,'  </SpecimenStartPoint>\r\n');

fprintf(fid,'  <SpecimenAngle>%1.0f</SpecimenAngle>\r\n',0); %specimen angle

fprintf(fid,'  <Comment>%s</Comment>\r\n',Commentstr); %comment

fprintf(fid,'  <Userfields>\r\n'); %userfields
for j=1:9 %write the Userfield entries
    fprintf(fid,'    <Userfield UserfieldID="Userfield %d">\r\n',j);
    fprintf(fid,'      <Value>\r\n      </Value>\r\n');
    fprintf(fid,'    </Userfield>\r\n');
end
fprintf(fid,'  </Userfields>\r\n'); %end userfields input

%write rows
for j=1:1
    fprintf(fid,'    <Row RowName="Row %d">\r\n',j);
    fprintf(fid,'      <KindOfMeasurement>Vickers</KindOfMeasurement>\r\n');
    fprintf(fid,'      <RowAngle>%1.0f</RowAngle>\r\n',0);
    fprintf(fid,'      <Status>\r\n      </Status>\r\n');
    fprintf(fid,'      <DateTime>%s</DateTime>\r\n',dt);
    fprintf(fid,'      <Method>%s</Method>\r\n',HVstr);
    fprintf(fid,'      <Objective>%s</Objective>\r\n',Objstr);
    fprintf(fid,'      <UseConversion>No</UseConversion>\r\n');
    fprintf(fid,'      <ConversionTable>\r\n      </ConversionTable>\r\n');
    fprintf(fid,'      <ConversionMaterial>\r\n      </ConversionMaterial>\r\n');
    fprintf(fid,'      <ConversionMethod>\r\n      </ConversionMethod>\r\n');
    fprintf(fid,'      <UseGeometryCorrection>No</UseGeometryCorrection>\r\n');
    fprintf(fid,'      <Shape>\r\n      </Shape>\r\n');
    fprintf(fid,'      <Curvature>\r\n      </Curvature>\r\n');
    fprintf(fid,'      <GeomCorrDiameter>\r\n      </GeomCorrDiameter>\r\n');
    fprintf(fid,'      <Angle>\r\n      </Angle>\r\n');
    fprintf(fid,'      <HardnessMin>0</HardnessMin>\r\n');
    fprintf(fid,'      <HardnessMax>0</HardnessMax>\r\n');
    fprintf(fid,'      <UseAutomaticIndentSpacing>No</UseAutomaticIndentSpacing>\r\n');
    fprintf(fid,'      <DistanceFromEdge>\r\n      </DistanceFromEdge>\r\n');
    fprintf(fid,'      <DistanceFactorAutomIndentSpacing>\r\n      </DistanceFactorAutomIndentSpacing>\r\n');
    fprintf(fid,'      <NumberOfIndents>\r\n      </NumberOfIndents>\r\n');
    fprintf(fid,'      <ZoomLevel>%1.0f</ZoomLevel>\r\n',zoom);
    fprintf(fid,'      <CircularLightUsed>No</CircularLightUsed>\r\n');
    fprintf(fid,'      <StartPoint>\r\n'); %use global start point
    fprintf(fid,'        <XAbs>%d</XAbs>\r\n',GlobalStartPoint(1));
    fprintf(fid,'        <YAbs>%d</YAbs>\r\n',GlobalStartPoint(2));
    fprintf(fid,'      </StartPoint>\r\n');
    
    %now do points
    for k=1:size(p,1)
        fprintf(fid,'      <Point PointID="%d">\r\n',k);
        fprintf(fid,'        <Hardness>\r\n        </Hardness>\r\n');
        fprintf(fid,'        <ImagePath>\r\n        </ImagePath>\r\n');
        fprintf(fid,'        <NPX>\r\n        </NPX>\r\n');
        fprintf(fid,'        <NPY>\r\n        </NPY>\r\n');
        fprintf(fid,'        <EPX>\r\n        </EPX>\r\n');
        fprintf(fid,'        <EPY>\r\n        </EPY>\r\n');
        fprintf(fid,'        <SPX>\r\n        </SPX>\r\n');
        fprintf(fid,'        <SPY>\r\n        </SPY>\r\n');
        fprintf(fid,'        <WPX>\r\n        </WPX>\r\n');
        fprintf(fid,'        <WPY>\r\n        </WPY>\r\n');
        fprintf(fid,'        <FocusPosition>\r\n        </FocusPosition>\r\n');
        fprintf(fid,'        <Diag>\r\n        </Diag>\r\n');
        fprintf(fid,'        <Diag1>\r\n        </Diag1>\r\n');
        fprintf(fid,'        <Diag2>\r\n        </Diag2>\r\n');
        fprintf(fid,'        <Classification>\r\n        </Classification>\r\n');
        fprintf(fid,'        <Status>\r\n        </Status>\r\n');
        fprintf(fid,'        <XAbs>%d</XAbs>\r\n',AbsMeasPnts(k,1));
        fprintf(fid,'        <YAbs>%d</YAbs>\r\n',AbsMeasPnts(k,2));
        fprintf(fid,'        <XRel>%0.3f</XRel>\r\n',RelMeasPnts(k,1));
        fprintf(fid,'        <YRel>%0.3f</YRel>\r\n',RelMeasPnts(k,2));
        fprintf(fid,'        <DateTime>\r\n        </DateTime>\r\n');
        fprintf(fid,'        <KindOfMeasurement>Vickers</KindOfMeasurement>\r\n');
        fprintf(fid,'        <Method>%s</Method>\r\n',HVstr);
        fprintf(fid,'        <Objective>%s</Objective>\r\n',Objstr);
        fprintf(fid,'        <ZoomLevel>%1.0f</ZoomLevel>\r\n',zoom);
        fprintf(fid,'        <UseConversion>No</UseConversion>\r\n');
        fprintf(fid,'        <ConversionTable>\r\n        </ConversionTable>\r\n');
        fprintf(fid,'        <ConversionMaterial>\r\n        </ConversionMaterial>\r\n');
        fprintf(fid,'        <ConversionMethod>\r\n        </ConversionMethod>\r\n');
        fprintf(fid,'        <ConversionValue>\r\n        </ConversionValue>\r\n');
        fprintf(fid,'        <UseGeometryCorrection>No</UseGeometryCorrection>\r\n');
        fprintf(fid,'        <Shape>\r\n        </Shape>\r\n');
        fprintf(fid,'        <Curvature>\r\n        </Curvature>\r\n');
        fprintf(fid,'        <GeomCorrDiameter>\r\n        </GeomCorrDiameter>\r\n');
        fprintf(fid,'        <Angle>\r\n        </Angle>\r\n');
        fprintf(fid,'        <User>%s</User>\r\n',User);
        fprintf(fid,'        <CircularLightUsed>No</CircularLightUsed>\r\n');
        for kk=1:3
            fprintf(fid,'        <AdditionalTestPointValue%d>\r\n        </AdditionalTestPointValue%d>\r\n',kk,kk);
        end
    
    
        fprintf(fid,'      </Point>\r\n');
    end
    fprintf(fid,'    </Row>\r\n');%end row
    
    
    
end
fprintf(fid,'</Specimen>'); %end specimen field
fclose(fid);


