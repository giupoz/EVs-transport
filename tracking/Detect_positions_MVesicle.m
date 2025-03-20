
% this code opens a movie in avi format and lets the user track a
% microvesicle or another object moving within the frame

close all
clear all

% ----- select and read avi movie -----

[FILENAME, PATHNAME] = uigetfile('*.avi', 'Load File .avi');
if (FILENAME == 0)
    disp('File not selected. Please select a valid file.');
    return;
end

filename = fullfile(PATHNAME, FILENAME);
data_s = VideoReader(filename); % movie data is read in variable data_s 
info=get(data_s) % info about the movie is displayed in command window
nFrames = data_s.NumberOfFrames % number of frames in the movie

Samp_1=rgb2gray(read(data_s,1)); % read first frame in variable Samp_1
h=figure(1); imshow(Samp_1); % display first frame in figure (1)
title('first frame'); colormap gray; axis square;
 
%--------------------------------------------------------
% SECTION for image properties
conv_fac=10; % conversion factor : 10 pixels makes 1 um
% bining=1; 
% magnefication=1;

%------------------- begins old version of LADAN code 
B = input('Introduce the number of the FIRST frame for analysis : ')
rep_1 = input('Do you want to process the movie to the end ? Y/N [Y]: ', 's')
if isempty(rep_1)
    rep_1 = 'Y';
end
if rep_1=='N'
    L= input('Introduce the number of the Last frame to be analized: ')
else 
    L= nFrames;    
end

Step= input('Introduce the number of STEP (skip) frames: ') % length of frame step

M0=rgb2gray(read(data_s,B));
h=figure(2); imshow(M0); hold on;
title(['Frame = ',num2str(B)]); colormap gray; axis square; % displays frame B in figure (2)

disp ('CLICK on MicroVesicle ')

[X_1 Y_1]= ginput(1)    

p_X(1)=floor(X_1);
p_Y(1)=floor(Y_1);
plot(p_X(1) ,p_Y(1),'xg'); 

SetFrame=B:Step:L;

hold off;

for k=2:numel(SetFrame);
    f=SetFrame(k);
    Mior=rgb2gray(read(data_s,f));

    h = figure(3) 
    imshow(Mior); hold on; 
    title(['Frame = ',num2str(f)]); colormap gray; axis square;
    plot(p_X(1) ,p_Y(1),'xg');
    hold on;
        disp ('CLICK on MicroVesicle ')
        [X Y]= ginput(1)    
        p_X(k)=floor(X); p_Y(k)=floor(Y);
    plot(p_X(k), p_Y(k),'or');
    hold off;
end

% plot the vsicle trajectory on the first image
Mior=rgb2gray(read(data_s,1));
h=figure(4); 
imshow(Mior); hold on
plot(p_X, flipud(p_Y),'.-r')
title(['Trajectory defined by = ',num2str(numel(SetFrame+1)), ' points (frames)']); colormap gray; axis square;
hold off

% save the coordinates of the trajectory nodes in file txt - for Excel
%P=[p_X', p_Y']
Name=[FILENAME,'MV.txt']
%mkdir([PATHNAME,'\analyzed'])
%save([PATHNAME,'\analyzed\',Name],'p_X', 'p_Y', '-ASCII', '-tabs');
% xlswrite([PATHNAME,'\analyzed\',Name],P)

P=p_X'
PP=p_Y'
%mkdir([PATHNAME,'\analyzed'])
fid=fopen([PATHNAME,'\analyzed12\',Name],'wt');
fprintf(fid,'%f\t%f\n',[P(:),PP(:)]');
fclose(fid); 


