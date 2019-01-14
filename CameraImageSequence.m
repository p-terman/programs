function status = CameraImageSequence(LED_on_filename,LED_on_path,dark,dark_path, data_folderName)

LED_location=[LED_on_path filesep LED_on_filename];
%LED=imread('/home/paul/Desktop/LED', 'jpeg');
LED=imread( LED_location , 'jpeg');

LED2=im2double(LED);
import matlab.io.*;
darkname= [dark_path filesep dark];
%darkname='/home/paul/Desktop/Yale_camera/dark/20141210-194225.fits';%specify for each set
dark=fitsread(darkname,'image');
%data_folderName='/home/paul/Desktop/Yale_camera/BD1/';
extension='*.fits';

concattedString=strcat(data_folderName, extension);
fileSet=dir(concattedString); 

for k = 1:length(fileSet);
    filename = strcat(data_folderName, fileSet(k).name);
    
   % fptr = fits.openFile(filename);
    S=filename(1:end-5);%remove .fits from name
    jFileName=strcat(S, '.jpeg');
    
    
    %dataStruct.(S) = fitsread(filename,'image'); %reads the fits file to the new name in the dataStruct
    data=fitsread(filename,'image');% this line is inplace of the one above. much faster not to save the dataStruct, and it causes crash with if many files
    %fits.closeFile(fptr);
    
    %sum=sum+dataStruct.(S);
    %sum=sum+data; % this is in place of the one above.
  %  fileSet(k).name

    cor=data-dark;
    pix_min;

    %tojpeg=enhanced;
    
    %if  max(tojpeg(:))>256;
     %   y=max(tojpeg(:))/256;
     %   tojpeg=enhanced/y;
    %end; %now every element of tojpeg should be less than 256 needed for jpeg saving
    
    m = flip(tojpeg ,1);
    
    I=double(zeros([size(tojpeg) 3]));
    I(:,:,1) = m+LED2;
    I(:,:,2) = LED2;
    I(:,:,3) = LED2;
    
    imshow(I)
    text(0,2000, fileSet(k).name(1:end-5) ,'Color',[1 1 0],'FontSize',20)
    hFrame = getframe(gca);
    imwrite(hFrame.cdata, jFileName, 'jpeg');
    
    %imwrite(I, jFileName);
end;

imageNames = dir(fullfile(data_folderName,'*.jpeg'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile(data_folderName,'timelapse.avi'));
outputVideo.FrameRate = 1;
open(outputVideo);

for ii = 1:length(imageNames);
   img = imread(fullfile(data_folderName,imageNames{ii}));
   writeVideo(outputVideo,img);
end;

close(outputVideo);

