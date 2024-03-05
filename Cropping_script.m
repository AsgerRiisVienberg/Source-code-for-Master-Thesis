clear all
close all
clc

currentFile = mfilename;

% Get the full path of the file
programPath = which(currentFile);
programFolder = fileparts(programPath);
addpath(programFolder);

% %let the user specify the path
dataPath = uigetdir;
cd(dataPath)



%overwrite settings
overwriteExsitingCroppedPhotos = false;

if overwriteExsitingCroppedPhotos
    disp('Overwriting exiting Cropped_Photos folder')
else
    disp('Overwriting disabled, will skip if the Cropped_Photos folder already exits')
end

%find the folder names
listing = dir;
folderNames = {listing([listing.isdir] & ~ismember({listing.name}, {'.', '..'})).name};

% display message
disp('The photos in the following folders will be cropped and saved')
for i = 1:numel(folderNames)
    fprintf('%s\n',folderNames{i});
end
fprintf('\n')

for i = 1:numel(folderNames)
    cd(dataPath)
    fprintf('\nCropping photos in %s \n',folderNames{i})
    currentFolder = folderNames{i};
    %save the current folder for the experiment data
    currentFolderPath = sprintf('%s\\%s',pwd,currentFolder);

    % Change the directory to the experiment danta
    cd(sprintf('%s',currentFolderPath))
 

    %Check if the photos folder is in the dir
    if ~exist('Photos', 'dir')
        disp('The /Photos folder does not exits... skipping')
        continue
    end

    % Creat the new folder for the cropped images
    newFolderName ='Cropped_Photos';
    if ~exist(newFolderName, 'dir')
        % If it doesn't exist, create it
        mkdir(newFolderName)
    elseif overwriteExsitingCroppedPhotos
        % If it already exists, display a message
        fprintf(['Folder "', newFolderName, '" already exists.']);
        fprintf('Overwriting...\n')
        rmdir(newFolderName,'s');
        mkdir(newFolderName);
    else
        disp('The Cropped_Photos folder already exits.')
        disp('checking if any photos are missing...')
    end

    %change to the folder with the uncropped images
    cd(sprintf('%s\\Photos',pwd))
    if exist('LocalMask.m','file') == 2
        mask = @(x) LocalMask(x);
    else
        mask = @(x) GreenMask(x);
    end


    % Define the file extension you want to search for
    fileExtension = '.jpg';

    % Use the 'dir' function to list all files in the current folder
    fileList = dir(fullfile(pwd, ['*' fileExtension]));

    % Extract the filenames and store them in a cell array
    fileNames = {fileList.name};
    numberOfPhotos = numel(fileNames);

    % Loop over all the file and crop them and creat the progress bar
    PhotosProcessed = 0;
    pb = CmdLineProgressBar('Prcessing photos...');
    for j = 1:numel(fileNames)
        currentFileName = fileNames{j};
        %check if the phot has already been processed
        cd(sprintf('%s\\%s',currentFolderPath,newFolderName))        
        if exist(sprintf('%s',currentFileName),"file") == 2
            PhotosProcessed = PhotosProcessed + 1;
            pb.print(PhotosProcessed,numberOfPhotos)
            continue
        end
        cd(sprintf('%s\\Photos',currentFolderPath))



        % Load the image and rotat it
        im = imread(currentFileName);
        im = imrotate(im,90);

        %filter for the red color
        imBinary = mask(im);
        %fill holes
        imBinary = imfill(imBinary,'holes');
        imBinary = imopen(imBinary,strel("square",3));
         % make sure that the this is a square
        imBinary = imdilate(imBinary,strel("square",80));
        % keep the smaller blobs that are closet to the biggest one
        imBinary = bwareafilt(imBinary,1);
        %find the bounding box
        stats = regionprops(imBinary,'BoundingBox');
        stats.BoundingBox(4) = stats.BoundingBox(4)*1.2;
        stats.BoundingBox(2) = stats.BoundingBox(2) - 0.1*stats.BoundingBox(4);



        %Check that only one box has been found
        if size(stats,1)~=1
            imcropped = imcrop(im);
            pause
        else
            %crop the image
            imcropped = imcrop(im,stats.BoundingBox);
        end

        %Change to the new folder and save the image then change back
        cd(sprintf('%s\\%s',currentFolderPath,newFolderName))
        imwrite(imcropped,currentFileName)
        cd(sprintf('%s\\Photos',currentFolderPath))

        % Update status of the process
        PhotosProcessed = PhotosProcessed + 1;
        pb.print(PhotosProcessed,numberOfPhotos)

    end

end

disp('Cropping complete')

