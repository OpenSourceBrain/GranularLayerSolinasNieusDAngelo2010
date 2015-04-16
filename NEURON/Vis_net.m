function Y = Vis_net(dd,Y)

if nargin == 1
    disp('Load data...')
    Y = VM_load(dd,'V');
end
% disp('Generate map...')
% netplot(dd);
% disp('Generate Vm plot of stimulated cells...')
% Vm_plot(Y);
% disp('Generate movie...')
% avi_file = [dd '/' dd '_movie'];
% net_movie(Y,avi_file);
% system(['ffmpeg -y -i ' avi_file '.avi -target pal-dvd -s 560x420 ' avi_file '.mp4.avi']);

Y = LFP_plot(Y,'VM');
