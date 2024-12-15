%tic
clc
clear

loops = 10;
x = linspace(0,2*pi, loops);
y1 = sin(x);

M(loops) = struct('cdata',[],'colormap',[]);


h = figure;
h.Visible = 'off';
axis manual
axis([0,2*pi, -1.5, 1.5])
hold on
for j = 1:loops
    scatter(x(j),y1(j),'filled',MarkerFaceColor = 'b')

    M(j) = getframe;
end
hold off

h.Visible = 'on';
movie(M);

% % Create a video file
% videoFile = VideoWriter('random_plot_movie', 'MPEG-4'); % Specify the file name and format
% open(videoFile); % Open the video file for writing
% 
% % Write each frame to the video file
% for k = 1:loops
%     writeVideo(videoFile, M(k));
% end
% 
% close(videoFile); % Close the video file

%toc

%%

%tic
pause(10);
% Generate data for the comet plot
t = linspace(0, 1, 4 );
y = sin(t); % Example datacomet

% Create a video file
videoFile = VideoWriter('comet_plot_movie', 'MPEG-4'); % Specify the file name and format
open(videoFile); % Open the video file for writing

f = figure;
hold on;
axis manual
axis([0, 200, -150, 50])

% Create the comet plot and capture frames
for k = 1:length(t)
    % Clear the current plot
    clf;

    % Plot the comet
    comet(u(1:50:12000));
    title('Comet Plot');
    xlabel('Time');
    ylabel('Value');

    % Capture the frame
    frame = getframe(gcf);
    writeVideo(videoFile, frame); % Write the frame to the video file
end
hold off

close(videoFile); % Close the video file
%toc