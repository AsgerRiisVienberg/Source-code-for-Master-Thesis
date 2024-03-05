clear all
close all
clc

%user settings

% %let the user specify the path
dataPath = 'C:\Users\Asger Riis Vienberg\OneDrive\DTU\MSc\Kanddiat Special\Eksperiments\Elastollan\Elastollan1190A_startup\Elstollan_Temp37_Strainrate0.01';
cd(dataPath)
% analyse phtos rutine
disp('Loading data')

% ---Read data from experiment----
%read the specimen information file to get the sample thickness
fileID = fopen(sprintf('%s\\%s',dataPath,'specimen_information.txt'),'r');
fgetl(fileID);
initialSampleThickness = str2double(fgetl(fileID));
fgetl(fileID);
initialSampleWidth = str2double(fgetl(fileID));
fgetl(fileID);
initialMarkerLength = str2double(fgetl(fileID));
fclose(fileID);

%Load the calibratin file
load("Photo_Calibration.mat")


% --- Load experimental data
dataTable = dataLoader(sprintf('%s\\%s',dataPath,'specimen.dat'));
% --- Load photo names and check the timing---

%change to the folder with the cropped images
cd(sprintf('%s\\Cropped_Photos',dataPath))
% find the photos in the folder
% Define the file extension to search for
fileExtension = '.jpg';
% Use the 'dir' function to list all files in the current folder
fileList = dir(fullfile(pwd, ['*' fileExtension]));
% Extract the filenames and store them in a cell array
fileNames = {fileList.name};
% Find the number of photos
numberOfPhotos = numel(fileNames);

%check that the first photo timestamp and the first data timestamp
% of day approx match
cd(sprintf('%s\\Photos',dataPath))
% Read image metadata
timeStampPhotos = zeros(size(fileNames));
for j = 1:numel(fileNames)
    timeTaken = datetime(imfinfo(fileNames{j}).DateTime,'Format', 'yyyy:MM:dd HH:mm:ss');
    timeStampPhotos(j) = hour(timeTaken) + minute(timeTaken)/60 + second(timeTaken)/3600;
end

%check the difference between the first part of the data and the photo
%in seconds
signalTimeDifference = abs(timeStampPhotos(1) - 24*dataTable.("Time Of Day [d]")(1))*3600;

if signalTimeDifference>60
    error('Timestamps between photo and data is more than 60 seconds apart!')
elseif signalTimeDifference>10
    warning('Timestamps between photo and data is more than 10 seconds apart!')
end

[photosTakenTime, photosTakenIndex,validPhotos,numberOfPulses] = PhotoTimes(dataTable,numberOfPhotos);
numberOfPhotos = numel(validPhotos);

%%
% check that timeings are uniform locate the timeing of the photos form the signal data


%check if the timestamps are even and if not ask the user for help
if numel(unique(round(diff(timeStampPhotos(validPhotos)),6)))>2
    plot(timeStampPhotos)
    x = input('Are the camara photo timestamps even? [y/n/pick]','s');
    close all
    if all(lower(x) == 'pick')
        validPhotos = PickPointsLinePlot(timeStampPhotos(validPhotos));
        photosTakenTime = photosTakenTime(validPhotos);
        photosTakenIndex = photosTakenIndex(validPhotos);
    elseif lower(x) == 'n'
        error('The timestamps are not even, please do something')
    end
end


%% ---analysins the strain in the photos---

%change back to the folder with the cropped images
cd(sprintf('%s\\Cropped_Photos',dataPath))

% Loop over all the and analyse them
deformation = zeros(numberOfPhotos,2);
sampleSliped = false;
PhotosProcessed = 0;
pb = CmdLineProgressBar('Processing photos...');
for i = 1:numberOfPhotos
    % Load the image
    currentFileName = fileNames{i};
    im = imread(currentFileName);

    % find the square and approxmiate its positiont
    if i == 1
        imAreaLast = ceil(sum(GreenMask(im),'all')*0.9);
    end
    % Load the image and rotat it


    ROI = GreenMask(im);
    % imshowpair(im,ROI)
    ROI = imfill(ROI,8,'holes');
    %save the original filter
    ROI_unfiltered = ROI;
    %filter all the smalle things out
    ROI = imopen(ROI,strel("square",20));
    % make sure that the this is a square
    ROI = imdilate(ROI,strel("square",10));

    % keep the largest blob
    ROI = bwareafilt(ROI,1);
    ROI_unfiltered = bwareaopen(ROI_unfiltered,1000,4);
    
    % imshowpair(im,ROI_unfiltered)
    imArea = sum(ROI(:));

    if 0.8*imAreaLast < imArea
        imAreaLast = imArea;
    else
        figure(1)
        cd(sprintf('%s\\Photos',dataPath))
        imshow(imread(currentFileName))
        cd(sprintf('%s\\Cropped_Photos',dataPath))
        disp('The Algorithem failed to find the square')
        x = input('Did the sample slide out of the jaws? [y/n]\n','s');
        if x == 'y' || x == 'Y'
            sampleSliped = true;
            deformation = deformation(deformation(:,1)~=0,:);
            photosTakenIndex = photosTakenIndex(1:(i-1));
            photosTakenTime = photosTakenTime(1:(i-1));            
            break
        else
            close all
            error('The alogritem cannot find the square')
        end
    end

    % find the max ferret properties
    ferret = bwferet(bwareafilt(ROI_unfiltered,1),'all');
    %rotat to the correct orientation

    if ferret.MinAngle>90
        rotAngle = ferret.MinAngle - 180;
    elseif ferret.MinAngle<-90
        rotAngle = ferret.MinAngle + 180;
    else
        rotAngle = ferret.MinAngle;
    end

    imsize = size(ROI);
    ROI = imrotate(ROI,rotAngle,'bicubic','crop');
    ROI_unfiltered = imrotate(ROI_unfiltered,rotAngle,'bicubic','crop');

    % keep the largest blob
    ROI = bwareafilt(ROI,1);
    ROI_unfiltered = bwareaopen(ROI_unfiltered,1000);
    
    stats = regionprops(ROI,'BoundingBox');
    bboxPoints = bbox2points(stats.BoundingBox);
    ROI = poly2mask(bboxPoints(:,1),bboxPoints(:,2),imsize(1),imsize(2));

    %find the edge and remove everything outside the ROI
    
    
    
    imedge = edge(imadjust(im2gray(im)),'canny',[0.01,0.05]);
    imedge = imrotate(imedge,rotAngle,'bicubic','crop');
    imedge(~ROI) = 0;
    imshowpair(imrotate(im,rotAngle,'bicubic','crop'),imedge)

    


    %messure elongation
    % messure using the green marking
    imtemp = zeros(imsize);
    index = ceil(stats.BoundingBox(1)+1.5/10*stats.BoundingBox(3)):...
            ceil(stats.BoundingBox(1)+8.5/10*stats.BoundingBox(3));
    imtemp(:,index) = ROI_unfiltered(:,index);
    sampleMessure = zeros(size(ROI_unfiltered,2),1);
    for j = 1:1:numel(sampleMessure)
        if any(imtemp(:,j) == 1)
            sampleMessure(j) = find(imtemp(:,j),1,'last') - find(imtemp(:,j),1,'first');
        end
    end
    sampleMessure = sampleMessure(sampleMessure>0);
    sampleMessure = sampleMessure(abs(sampleMessure-mean(sampleMessure))<2*std(sampleMessure));
    deformation(i,1) = mean(sampleMessure(any(sampleMessure > 0,2)));

    % messure contration
    imtemp = zeros(size(imedge));
    index = ceil(stats.BoundingBox(2)+1/10*stats.BoundingBox(4)):...
            ceil(stats.BoundingBox(2)+9/10*stats.BoundingBox(4));
    imtemp(index,:) = imedge(index,:);
    
    sampleMessure = zeros(size(imedge,1),1);
    for j = 1:1:numel(sampleMessure)
        if any(imedge(j,:) == 1)
            sampleMessure(j) = find(imedge(j,:),1,'last') - find(imedge(j,:),1,'first');
        end
    end
    %filter
    sampleMessure = sampleMessure(sampleMessure > 0);
    sampleMessure = sampleMessure( abs(sampleMessure-mean(sampleMessure))<1*std(sampleMessure) );
    deformation(i,2) = mean(sampleMessure(any(sampleMessure > 0,2)));

    PhotosProcessed = PhotosProcessed + 1;
    pb.print(PhotosProcessed,numberOfPhotos)
end



%% Save data
%convert displacement from pixsels to mm
close all
displacements(:,2) = deformation(:,2)*photoCalibration
plot(displacements)
%%
close all
point = 278
displacements(point,1) = (displacements(point-1,1)+displacements(point+1,1))/2;
plot(displacements())
%%


%%
% dataPath = 'C:\Users\Asger Riis Vienberg\OneDrive\DTU\MSc\Kanddiat Special\Eksperiments\Elastollan\Pre-test 1\Elastollan_Temp37_Strainrate0.01';
% dataFilePath = sprintf('%s\\%s',dataPath,'specimen.dat');
% data1= dataLoader(dataFilePath);
% 
% 
% dataPath = 'C:\Users\Asger Riis Vienberg\OneDrive\DTU\MSc\Kanddiat Special\Eksperiments\Elastollan\Pre-test 1\Elastollan_Temp37_Strainrate0.1';
% dataFilePath = sprintf('%s\\%s',dataPath,'specimen.dat');
% data2= dataLoader(dataFilePath);

%%

plot(data1.("Ch 1 Displacement [mm]"),data1.("Ch 1 Load 4B [N]"))
hold on
plot(data2.("Ch 1 Displacement [mm]"),data2.("Ch 1 Load 4B [N]"))

%%
close all
plot(dataTable.("Ch 1 Load 4B [N]")(end-200000:end))
%%
close all
zeroforce = mean(dataTable.("Ch 1 Load 4B [N]")(end-100000:end))


dataTable.("Ch 1 Load 4B [N]") = dataTable.("Ch 1 Load 4B [N]")-zeroforce;
plot(dataTable.("Ch 1 Load 4B [N]"))
yline(0)


%%

% photosTakenIndex(69:end) = photosTakenIndex(69:end)+2048;

close all
plot(dataTable.("Time [s]")(photosTakenIndex),(displacements(:,1)-displacements(1,1))/(max(displacements(:,1))-displacements(1,1))*0.75,'.-')
hold on
plot(dataTable.("Time [s]"),(dataTable.("Ch 1 Displacement [mm]")+120)/max(dataTable.("Ch 1 Displacement [mm]")+120),'.-')

%%
close all
plot(displacements(:,1),dataTable.("Ch 1 Load 4B [N]")(photosTakenIndex),'.-')
%%
f = fit(brushedData(:,1),brushedData(:,2),'poly1')
plot(f,brushedData(:,1),brushedData(:,2))
initialMarkerLength = -f.p2/f.p1

% f = fit(brushedData(:,1),brushedData(:,2),'poly2')
% plot(f,brushedData(:,1),brushedData(:,2))
% initialMarkerLength = (-f.p2+sqrt(f.p2^2-4*f.p1*f.p3))/(2*f.p1)
% (-f.p2-sqrt(f.p2^2-4*f.p1*f.p3))/(2*f.p1)


%%
close all
plot(displacements(:,2),dataTable.("Ch 1 Load 4B [N]")(photosTakenIndex),'.-')

%%
% f = fit(brushedData(:,1),brushedData(:,2),'poly1')
% plot(f,brushedData(:,1),brushedData(:,2))
% intersect = -f.p2/f.p1
% initialSampleWidth


f = fit(brushedData(:,1),brushedData(:,2),'poly2')
plot(f,brushedData(:,1),brushedData(:,2))
(-f.p2+sqrt(f.p2^2-4*f.p1*f.p3))/(2*f.p1)
intersect = (-f.p2-sqrt(f.p2^2-4*f.p1*f.p3))/(2*f.p1)
initialSampleWidth

%%
displacements(displacements(:,2)~=initialSampleWidth,2)=displacements(displacements(:,2)~=initialSampleWidth,2) - (intersect-initialSampleWidth)

%%

plot(displacements(:,2),dataTable.("Ch 1 Load 4B [N]")(photosTakenIndex),'.-')
x = linspace(5,initialSampleWidth);
hold on
plot(x,f(x))

%%

%%
strains = log(displacements./[initialMarkerLength initialSampleWidth])
crossSection = displacements(:,2) .* (exp(strains(:,2))*initialSampleThickness);
stress = dataTable.("Ch 1 Load 4B [N]")(photosTakenIndex)./crossSection;
figure()
plot(strains(:,1),stress(:),'.-')
xline(0)
yline(0)
%%
sampleData = table(displacements(:,1),displacements(:,2), strains(:,1),strains(:,2),crossSection(:),stress(:), ...
    'VariableNames', {'Vertical displacement [mm]','Horizontal displacement [mm]','Strain-11', 'Strain-22', 'Cross section [mm2]','Stress [MPa]'});

%change the directory
cd(dataPath);

%combing the tables at the photos taken and save the table
data = [dataTable(photosTakenIndex(:),:),sampleData];

%%
%interpolate the strain values between the photos
% Interpolated the messured data

%%
idxInterp = photosTakenIndex(1):photosTakenIndex(end);
displacementsInterp = zeros(size(idxInterp,2),2);
displacementsInterp(:,1) = interp1(dataTable.("Time [s]")(photosTakenIndex),displacements(:,1),dataTable.("Time [s]")(idxInterp),'linear');
displacementsInterp(:,2) = interp1(dataTable.("Time [s]")(photosTakenIndex),displacements(:,2),dataTable.("Time [s]")(idxInterp),'linear');

%%
close all


plot(displacementsInterp(:,1),'.-')

%%
close all
plot(dataTable.("Time [s]")(photosTakenIndex),displacements(:,1),'r.')
hold on
plot(dataTable.("Time [s]")(idxInterp),displacementsInterp(:,1),'-')
plot(dataTable.("Time [s]"),dataTable.("Ch 1 Load 4B [N]"))
plot(dataTable.("Time [s]")(photosTakenIndex),dataTable.("Ch 1 Load 4B [N]")(photosTakenIndex),'b.')
%%
plot(displacementsInterp(:,1),dataTable.("Ch 1 Load 4B [N]")(idxInterp))
hold on
plot(displacements(:,1),dataTable.("Ch 1 Load 4B [N]")(photosTakenIndex),'r.')

%%
close all
plot(dataTable.("Time [s]")(idxInterp),dataTable.("Ch 1 Load 4B [N]")(idxInterp))
hold on
plot(dataTable.("Time [s]")(photosTakenIndex),dataTable.("Ch 1 Load 4B [N]")(photosTakenIndex),'r.')
%%
plot(displacementsInterp)
hold on
plot(dataTable.("Ch 1 Load 4B [N]")(idxInterp))
%%
strainsInterp = log(displacementsInterp./[initialMarkerLength initialSampleWidth]);
crossSectionIntrep = displacementsInterp(:,2) .* (exp(strainsInterp(:,2))*initialSampleThickness);
stressInterp = dataTable.("Ch 1 Load 4B [N]")(idxInterp)./crossSectionIntrep(:);
close all
plot(strainsInterp(:,1),stressInterp)
%%
%covent back to table
sampleDataInterp = table(displacementsInterp(:,1),displacementsInterp(:,2), strainsInterp(:,1),strainsInterp(:,2),crossSectionIntrep,stressInterp, ...
    'VariableNames', {'Vertical displacement [mm]','Horizontal displacement [mm]','Strain-11', 'Strain-22', 'Cross section [mm2]','Stress [MPa]'});
dataInterp = [dataTable(idxInterp,:),sampleDataInterp];



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataTable = dataLoader(file_path, dataLines)
    %IMPORTFILE Import data from a text file
    %  DataLoader read the datafile from Structural lab.
    %  reads data from text file FILE_PATH for the default selection.
    %  Returns the data as a table.
    %  Table Colum ames are defined form line 3 of the inputfile unless
    % otherwise specified.
    
    %check that the path is absolut
    if ~startsWith(file_path, 'C:\')
        error('Please provide an absolut file for the folder or the specific file')
    end

    %Check whetear the file or folder is specified
    if ~endsWith(file_path,'\specimen.dat')
        file_path = append(file_path,'\specimen.dat');
    end
    
    % find names of varibles    
    if nargin < 2
        dataLines = 3;
    end

    fileID = fopen(file_path, 'r');
    
    % Skip the first 7 lines (header)
    for i = 1:dataLines
        fgetl(fileID);
    end
    
    % Read the variable names from line 8
    headerLine = fgetl(fileID);
    variableDescription = string(strsplit(headerLine, '\t')); % Assuming tab delimiter, adjust if needed
    unitLine = fgetl(fileID);
    variableUnits = string(strsplit(unitLine, '\t'));
    
    
    % Combine the strings
    variableNames = cellfun(@(variableDescription, variableUnits) [variableDescription ' [' variableUnits ']'], variableDescription, variableUnits, 'UniformOutput', false);
    
    numberOfVariables = length(variableNames); % This should dynamically determine the number of variables
    % Close the file
    fclose(fileID);
    
    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", numberOfVariables);
    
    % Specify range and delimiter
    opts.DataLines = [dataLines+1, Inf];
    opts.Delimiter = ["", "\t"];
    
    % Specify column names and types
    opts.VariableNames = variableNames;
    opts.PreserveVariableNames = true;
    opts.VariableTypes = repmat({'double'},1,numberOfVariables);
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts,"DecimalSeparator", ",");
    opts = setvaropts(opts,"ThousandsSeparator", ".");
    
    % Import the data
    dataTable = readtable(file_path, opts);
    nanRows = any(ismissing(dataTable),2);
    dataTable(nanRows, :) = [];
    
    % Clear temporary variables
    clear opts
end

function [peakCenterTime, peakCenterIndex,validPhotos,numberOfPulses] = PhotoTimes(dataTable,numberOfPhotos)
binarySignal = dataTable.("BNC 4A [unitless]");

binarySignal(1:10)=1;
% Find the indices of the rising and falling edges of the peaks
risingEdges = find(diff(binarySignal) == 1);
fallingEdges = find(diff(binarySignal)== -1)+1;

numberOfPulses = length(fallingEdges);

%calculate the peak width and filte the ones to small to trigger the camera
peakWidth = risingEdges - fallingEdges;

risingEdges = risingEdges(peakWidth >= 4);
fallingEdges = fallingEdges(peakWidth >= 4);

% Find the center of each peak and calculate peak width
validPhotos = 1:min(numberOfPhotos,numel(fallingEdges));
peakCenterTime = zeros(size(validPhotos));
peakCenterIndex = zeros(size(validPhotos));
for i = 1:min(numberOfPhotos,numel(fallingEdges))
    idx = ceil(median(fallingEdges(i):risingEdges(i)));
    peakCenterTime(i) = dataTable.("Time [s]")(idx);
    peakCenterIndex(i) = idx;
end

end


function validIndices = PickPointsLinePlot(array)
    % Plot the vector as a line plot
    figure;
    plot(array, 'o-');
    title('Array Line Plot');
    xlabel('Index');
    ylabel('Value');

    % Allow the user to interactively pick two points on the plot
    disp('Click on two points on the plot.');
    pickIndices = zeros(2,1);
    for i = 1:2
        % Wait for user input
        [x, ~] = ginput(1);

        % Round the index to the nearest integer
        index = round(x);

        % Plot a red circle around the selected point
        hold on;
        scatter(index, array(index), 'r', 'filled');
        hold off;
        
        % save the picked points
        pickIndices(i) = index;
    end
    pause(2)
    close all
    validIndices = [1:min(pickIndices) max(pickIndices):length(array)];
end


function [BW,maskedRGBImage] = GreenMask(RGB)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder app. The colorspace and
%  range for each channel of the colorspace were set within the app. The
%  segmentation mask is returned in BW, and a composite of the mask and
%  original RGB images is returned in maskedRGBImage.

% Auto-generated by colorThresholder app on 27-Dec-2023
%------------------------------------------------------


% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.200;
channel1Max = 0.563;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.048;
channel2Max = 0.510;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.294;
channel3Max = 0.992;

% Create mask based on chosen histogram thresholds
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end


