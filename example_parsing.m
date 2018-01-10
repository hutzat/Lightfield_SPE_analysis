%% creating all the spectral series objects. 


dot='0_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_0_degree=parse_LF_data('0_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_0_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_0_degree.contour_spectra(ylim_spec)


dot='30_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_30_degree=parse_LF_data('30_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_30_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_30_degree.contour_spectra(ylim_spec)


dot='60_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_60_degree=parse_LF_data('60_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_60_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_60_degree.contour_spectra(ylim_spec)


dot='90_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_90_degree=parse_LF_data('90_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_90_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_90_degree.contour_spectra(ylim_spec)


dot='120_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_120_degree=parse_LF_data('120_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_120_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_120_degree.contour_spectra(ylim_spec)

dot='150_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_150_degree=parse_LF_data('150_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_150_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_150_degree.contour_spectra(ylim_spec)

dot='180_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_180_degree=parse_LF_data('180_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_180_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_180_degree.contour_spectra(ylim_spec)


dot='210_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_210_degree=parse_LF_data('210_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_210_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_210_degree.contour_spectra(ylim_spec)


dot='240_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_240_degree=parse_LF_data('240_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_240_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_240_degree.contour_spectra(ylim_spec)



dot='270_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_270_degree=parse_LF_data('270_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_270_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_270_degree.contour_spectra(ylim_spec)



dot='300_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_300_degree=parse_LF_data('300_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_300_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_300_degree.contour_spectra(ylim_spec)



dot='330_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_330_degree=parse_LF_data('330_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_330_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_330_degree.contour_spectra(ylim_spec)



dot='360_degree';

%% reading in data,creating video, and creating contour plot.
ylim_spec=[2330,2365];
DotB_360_degree=parse_LF_data('360_degree.spe','dummy.ht2','dummy.ht3',1E6,'asd')
DotB_360_degree.create_spectral_series_video(strcat(dot,'.avi'))
DotB_360_degree.contour_spectra(ylim_spec)


%% create plot of the spectra as a a function of polarization degree.
P=addprop(DotB_0_degree,'avg_frames')
P=addprop(DotB_30_degree,'avg_frames')
P=addprop(DotB_60_degree,'avg_frames')
P=addprop(DotB_90_degree,'avg_frames')
P=addprop(DotB_120_degree,'avg_frames')
P=addprop(DotB_150_degree,'avg_frames')
P=addprop(DotB_180_degree,'avg_frames')
P=addprop(DotB_210_degree,'avg_frames')
P=addprop(DotB_240_degree,'avg_frames')
P=addprop(DotB_270_degree,'avg_frames')
P=addprop(DotB_300_degree,'avg_frames')
P=addprop(DotB_330_degree,'avg_frames')
P=addprop(DotB_360_degree,'avg_frames')

DotB_0_degree.avg_frames=mean(DotB_0_degree.frames(:,48:50),2);
DotB_30_degree.avg_frames=mean(DotB_30_degree.frames(:,3:5),2);
DotB_60_degree.avg_frames=mean(DotB_60_degree.frames(:,5:8),2);
DotB_90_degree.avg_frames=mean(DotB_90_degree.frames(:,10:12),2);
DotB_120_degree.avg_frames=mean(DotB_120_degree.frames(:,20:25),2);
DotB_150_degree.avg_frames=mean(DotB_150_degree.frames(:,20:25),2);
DotB_180_degree.avg_frames=mean(DotB_180_degree.frames(:,40:45),2);
DotB_210_degree.avg_frames=mean(DotB_210_degree.frames(:,1:10),2);
DotB_240_degree.avg_frames=mean(DotB_240_degree.frames(:,40:45),2);
DotB_270_degree.avg_frames=mean(DotB_270_degree.frames(:,10:15),2);
DotB_300_degree.avg_frames=mean(DotB_300_degree.frames(:,20:40),2);
DotB_330_degree.avg_frames=mean(DotB_330_degree.frames(:,40:50),2);
DotB_360_degree.avg_frames=mean(DotB_360_degree.frames(:,10:20),2);

%% fitting the polarization dependent spectra with the sum of two Lorentzians.


limx=[2340,2355];

figure('pos',[10,10,900,900])
subplot(3,5,1)
plot(DotB_0_degree.calib_eV,DotB_0_degree.avg_frames)
xlim(limx)
title('\Theta = 0')

subplot(3,5,2)
plot(DotB_30_degree.calib_eV,DotB_30_degree.avg_frames)
xlim(limx)
title('\Theta = 30')


subplot(3,5,3)
plot(DotB_60_degree.calib_eV,DotB_60_degree.avg_frames)
xlim(limx)
title('\Theta = 60')


subplot(3,5,4)
plot(DotB_90_degree.calib_eV,DotB_90_degree.avg_frames)
xlim(limx)
title('\Theta = 90')


subplot(3,5,5)
plot(DotB_120_degree.calib_eV,DotB_120_degree.avg_frames)
xlim(limx)
title('\Theta = 120')


subplot(3,5,6)
plot(DotB_150_degree.calib_eV,DotB_150_degree.avg_frames)
xlim(limx)
title('\Theta = 150')


subplot(3,5,7)
plot(DotB_180_degree.calib_eV,DotB_180_degree.avg_frames)
xlim(limx)
title('\Theta = 180')


subplot(3,5,8)
plot(DotB_210_degree.calib_eV,DotB_210_degree.avg_frames)
xlim(limx)
title('\Theta = 210')


subplot(3,5,9)
plot(DotB_240_degree.calib_eV,DotB_240_degree.avg_frames)
xlim(limx)
title('\Theta = 240')


subplot(3,5,10)
plot(DotB_270_degree.calib_eV,DotB_270_degree.avg_frames)
xlim(limx)
title('\Theta = 270')


subplot(3,5,11)
plot(DotB_300_degree.calib_eV,DotB_300_degree.avg_frames)
xlim(limx)
title('\Theta = 300')



subplot(3,5,12)
plot(DotB_330_degree.calib_eV,DotB_330_degree.avg_frames)
xlim(limx)
title('\Theta = 330')


subplot(3,5,13)
plot(DotB_360_degree.calib_eV,DotB_360_degree.avg_frames)
xlim(limx)
title('\Theta = 360')



%% creating a polar plot of the average intensity per frame over the whole 
%% time series. 
% I don't think there is any fine-structre state in DotG. 
P=addprop(DotB_0_degree,'avg_int')
P=addprop(DotB_30_degree,'avg_int')
P=addprop(DotB_60_degree,'avg_int')
P=addprop(DotB_90_degree,'avg_int')
P=addprop(DotB_120_degree,'avg_int')
P=addprop(DotB_150_degree,'avg_int')
P=addprop(DotB_180_degree,'avg_int')
P=addprop(DotB_210_degree,'avg_int')
P=addprop(DotB_240_degree,'avg_int')
P=addprop(DotB_270_degree,'avg_int')
P=addprop(DotB_300_degree,'avg_int')
P=addprop(DotB_330_degree,'avg_int')
P=addprop(DotB_360_degree,'avg_int')

DotB_0_degree.avg_int=sum(sum(DotB_0_degree.frames(290:450,:)))/size(DotB_0_degree.frames(290:450,:),2); 
DotB_30_degree.avg_int=sum(sum(DotB_30_degree.frames(290:450,:)))/size(DotB_30_degree.frames(290:450,:),2);
DotB_60_degree.avg_int=sum(sum(DotB_60_degree.frames(290:450,:)))/size(DotB_60_degree.frames(290:450,:),2);
DotB_90_degree.avg_int=sum(sum(DotB_90_degree.frames(290:450,:)))/size(DotB_90_degree.frames(290:450,:),2);
DotB_120_degree.avg_int=sum(sum(DotB_120_degree.frames(290:450,:)))/size(DotB_120_degree.frames(290:450,:),2);
DotB_150_degree.avg_int=sum(sum(DotB_150_degree.frames(290:450,:)))/size(DotB_150_degree.frames(290:450,:),2);
DotB_180_degree.avg_int=sum(sum(DotB_180_degree.frames(290:450,:)))/size(DotB_180_degree.frames(290:450,:),2);
DotB_210_degree.avg_int=sum(sum(DotB_210_degree.frames(290:450,:)))/size(DotB_210_degree.frames(290:450,:),2);
DotB_240_degree.avg_int=sum(sum(DotB_240_degree.frames(290:450,:)))/size(DotB_240_degree.frames(290:450,:),2);
DotB_270_degree.avg_int=sum(sum(DotB_270_degree.frames(290:450,:)))/size(DotB_270_degree.frames(290:450,:),2);
DotB_300_degree.avg_int=sum(sum(DotB_300_degree.frames(290:450,:)))/size(DotB_300_degree.frames(290:450,:),2);
DotB_330_degree.avg_int=sum(sum(DotB_330_degree.frames(290:450,:)))/size(DotB_330_degree.frames(290:450,:),2);
DotB_360_degree.avg_int=sum(sum(DotB_360_degree.frames(290:450,:)))/size(DotB_360_degree.frames(290:450,:),2);

theta=[0,30,60,90,120,150,180,210,240,270,300,330,360]/360*2*pi;
int=[DotB_0_degree.avg_int,DotB_30_degree.avg_int,DotB_60_degree.avg_int,DotB_90_degree.avg_int,DotB_120_degree.avg_int,DotB_150_degree.avg_int,DotB_180_degree.avg_int,DotB_210_degree.avg_int,DotB_240_degree.avg_int,DotB_270_degree.avg_int,DotB_300_degree.avg_int,DotB_330_degree.avg_int,DotB_360_degree.avg_int]
figure()
polar(theta,int)
title('Polarization Dependent Intensity Dot B')
set(gca,'fontsize',14)




