
load('HH_20170206135645_4.mat');
% load('NoiseFile.mat') ; 

Data_filtered = zeros(size(Data_out));
Data_filtered(3:end,:) = Data_out(3:end,:) - 2*Data_out(2:end-1,:) + Data_out(1:end-2,:);

% old output only
%Data_filtered = Data_out;



% img_file="C:\Users\14765\Documents\MATLAB\TUDelft\RadarI\practicum1\alldata";
%  time_ind = -500:1:500;
%  hfig=figure;
% 
%  imagesc(time_ind,range,db(abs(Data_filtered')))
%  colorbar
%  colormap('jet');
%  set(gca,'ydir','norm')
%  set(gca,'clim',[10,100]) 
% 
%  xlabel('Slow time, ms')
%  ylabel('Range, m')
%  title_str = "Initial Image Display";
%  title(['{',title_str,'}'])
%  print(hfig,'-dpng',img_file);


% N_Doppler=512; j=2;
% start_time=1+N_Doppler*(j-1);
% x=Data_filtered(start_time:start_time+N_Doppler-1,:);
% RD=fftshift(fft(x, N_Doppler),1);
% frequency=[-500:1000/(N_Doppler+1):500]; % how this has to be changed for diff PRF?
% 
% hfig=figure;
% imagesc(frequency,range,db(abs(RD')))
% colorbar
% colormap('jet');
% set(gca,'ydir','norm')
% set(gca,'clim',[10,140]) % If you do not see the range-Doppler plane similar to slide 10,
%  % comment (or edit) the codeline set(gca,'clim',[10,70])
% xlabel('Doppler frequency, Hz')
% ylabel('Range, m')
% title_str = "Range-Doppler map";
% title(['{',title_str,' 1ms, burst ',num2str(j),'}'])

PRI = 1;
PRI_0 = 1e-3;
PRF = 1e3;
fc = 3.315e9; 
c = 3e8;
lambda = c / fc;

imgDir_video1 = 'C:\Users\14765\Documents\MATLAB\TUDelft\RadarI\LAB1\';
N_Doppler=512;

frequency=-PRF/(2*PRI):(PRF/PRI)/(N_Doppler+1):PRF/(2*PRI); % how this has to be changed for diff PRF?
velocity = (lambda / 2) * frequency;

name = ['range_doppler_', num2str(N_Doppler), 'Ndopper_', num2str(PRI), 'PRI_filtered'];

video_file=[imgDir_video1,name,'.avi'];
writerObj = VideoWriter(video_file);
open(writerObj);

for j = 1:59
    start_time=1+N_Doppler*(j-1);
    if (start_time+PRI*N_Doppler-1 > 30720) 
        break; 
    end
    x=Data_filtered(start_time:PRI:start_time+PRI*N_Doppler-1,:);
    RD=fftshift(fft(x, N_Doppler),1);
    
    hfig = figure('Visible','off');
    % imagesc(frequency,range,db(abs(RD')));
    imagesc(velocity, range, db(abs(RD')));
    colorbar
    set(gca,'ydir','norm')
    set(gca,'clim',[10,150]) % If you do not see the range-Doppler plane similar to slide 10,
    colormap("jet")
    % xlabel('Doppler frequency, Hz')
    xlabel('Radial velocity, m/s')

    ylabel('Range, m')
    title_str = 'Range Doppler';
    title(['{',title_str,' 1ms, burst ',num2str(j),'}'])
    drawnow;
    frame = getframe(hfig);
    writeVideo(writerObj,frame);
    close(hfig)
end
close(writerObj);
