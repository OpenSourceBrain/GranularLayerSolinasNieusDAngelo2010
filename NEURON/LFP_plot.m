function Y = LFP_plot(Y,VM,GABA,LFP,LFP_G,mode_all)

% % Matrices of the coordinates in top view
% xmin = min(Y.goc_xyz(:,1));
% xmax = max(Y.goc_xyz(:,1));
% ymin = min(Y.goc_xyz(:,2));
% ymax = max(Y.goc_xyz(:,2));
% dxy = 1e-6;
% [XI,YI] = meshgrid(xmin:dxy:xmax, ymin:dxy:ymax);
% 
% Y.dt = 0.025; % ms
% Y.tstop = max(max(max(Y.goc_vm(:,1,:),[],3)),max(max(Y.goc_vm(:,1,:),[],3)));
% Y.t = 0:Y.dt:tstop;
% if nargin == 1
%     % LFP matrix 3Dim: 1 time or vm value, 2 time or vm, 3 cell
%     for igoc = 1:size(Y.goc_vm,3)
%         notnan = ~isnan(Y.goc_vm(:,1,igoc));
%         VM(:,igoc) = interp1(Y.goc_vm(notnan,1,igoc),Y.goc_vm(notnan,2,igoc),t);
%     end
% 
%     for it = find(t>50 & t<100) %length(t)
%         if mod(it,100)==0
%             disp(it)
%         end
%         LFP(:,:,it) = griddata(Y.goc_xyz(:,1),Y.goc_xyz(:,2),VM(it,:),XI,YI);
%     end
% end
% 
% % Open avi file
% mov = avifile([Y.dd '/' Y.dd '.avi']);
% 
% 
% for it = 2000:size(LFP,3)
%     surf(XI,YI,LFP(:,:,it))
% %     axis equal
%     shading interp
%     zlim([-70 -20])
%     caxis([-70 -20])
%     title(['Time ' sprintf('%05.3f',t(it)+50) ' ms'])
%     disp(['Time ' sprintf('%05.3f',t(it)+50) ' ms'])
%     view(2)
%     pause(0.01)
% %     print('-djpeg90',[Y.dd '/' Y.dd '.jpg'])
%     F = getframe(gcf);
%     mov = addframe(mov,F);    
% end
% mov = close(mov);


% Matrices of the coordinates in top view
xmin = min(Y.grc_xyz(:,1));
xmax = max(Y.grc_xyz(:,1));
ymin = min(Y.grc_xyz(:,2));
ymax = max(Y.grc_xyz(:,2));
dxy = 1e-6;
[XI,YI] = meshgrid(xmin:dxy:xmax, ymin:dxy:ymax);

Y.dt = 0.1; % ms
Y.tstop = max(max(max(Y.grc_vm(:,1,:),[],3)),max(max(Y.grc_vm(:,1,:),[],3)));
vm_min = min(min(min(Y.grc_vm(:,2,:),[],3)),min(min(Y.grc_vm(:,2,:),[],3)));
Y.t = 200:Y.dt:Y.tstop-100;
%Y.t = 0:Y.dt:200;
length(Y.t)
if nargin == 1
    for igrc = 1:size(Y.grc_vm,3)
        if mod(igrc,100)==0
            disp(igrc)
        end
        %igrc
        notnan = ~isnan(Y.grc_vm(:,1,igrc));
        vm = Y.grc_vm(notnan,:,igrc);
        notdup = ([diff(vm(:,1)); 1] > 0);
        %size(notdup)
        vm = vm(notdup,:);
        VM(:,igrc) = interp1(vm(:,1),vm(:,2),Y.t);

        notnan = ~isnan(Y.grc_gaba(:,1,igrc));
        vm = Y.grc_gaba(notnan,:,igrc);
        notdup = ([diff(vm(:,1)); 1] > 0);
        %size(notdup)
        vm = vm(notdup,:);
        GABA(:,igrc) = interp1(vm(:,1),vm(:,2),Y.t);
    end
    LFP = [];
    LFP_G = [];
end
 
if nargin == 3
    for k = 1 %:10
        %grc_idx(k).d = find(Y.grc_xyz(:,3)>-50e-6+(k-1)*10e-6 &Y.grc_xyz(:,3)<-50e-6+(k)*10e-6);
        grc_idx(k).d = find(Y.grc_xyz(:,3)>-20e-6 &Y.grc_xyz(:,3)<20e-6)
        %    size(grc_idx(k).d)
    end   % LFP matrix 3Dim: 1 time or vm value, 2 time or vm, 3 cell

    for it = 1:length(Y.t) %find(t>50 & t<100)
        if mod(it,100)==0
            disp(it)
        end
        %        for k = 1 %:10
        LFP(:,:,it) = griddata(Y.grc_xyz(grc_idx(1).d,1),Y.grc_xyz(grc_idx(1).d,2),VM(it,grc_idx(1).d),XI,YI);
        LFP_G(:,:,it) = griddata(Y.grc_xyz(grc_idx(1).d,1),Y.grc_xyz(grc_idx(1).d,2),GABA(it,grc_idx(1).d),XI,YI);
    end
    %    end
end

if nargin == 5
    % Open avi file
    mov = avifile([Y.dd '/' Y.dd '.avi']);
    
    figure(1)
    for it = 1:size(LFP,3)
        if mod(it,1000)==0
            disp(Y.t(it))
        end
        pcolor(LFP(:,:,it))
        shading interp
        lighting phong        
    %    zlim([vm_min -20])
        caxis([-75 -20])
        colorbar
        title(['Time ' sprintf('%05.3f',Y.t(it)) ' ms'])
        disp(['Time ' sprintf('%05.3f',Y.t(it)) ' ms'])
    %    view(2)
    %hold off
    %refresh(gcf)
    %drawnow update
    pause(.01)
    %close
    %     print('-djpeg90',[Y.dd '/' Y.dd '.jpg'])
       % F = getframe(gcf);
       % mov = addframe(mov,F);
    end
    mov = close(mov);
    
%     % Open avi file
%     mov = avifile([Y.dd '/' Y.dd 'GABA.avi']);
%     
%     figure(2)
%     for it = 1:size(LFP_G,3)
%         if mod(it,1000)==0
%             disp(Y.t(it))
%         end
%         pcolor(LFP_G(:,:,it))
%         shading interp
%     %    zlim([vm_min -20])
%         caxis([400 1600])
%         colorbar
%         title(['Time ' sprintf('%05.3f',Y.t(it)) ' ms'])
%         disp(['Time ' sprintf('%05.3f',Y.t(it)) ' ms'])
%     %    view(2)
%         pause(0.01)
%     %     print('-djpeg90',[Y.dd '/' Y.dd '.jpg'])
%         F = getframe(gcf);
%         mov = addframe(mov,F);
%     end
%     mov = close(mov);
end

if nargin == 2
    % VM & GABA
    VM = [];
    if ~isfield(Y,'VM')
        for igrc = 1:size(Y.grc_vm,3)
            if mod(igrc,100)==0
                disp(igrc)
            end
            %igrc
            notnan = ~isnan(Y.grc_vm(:,1,igrc));
            vm = Y.grc_vm(notnan,:,igrc);
            notdup = ([diff(vm(:,1)); 1] > 0);
            %size(notdup)
            vm = vm(notdup,:);
            Y.VM(:,igrc) = interp1(vm(:,1),vm(:,2),Y.t);

            if isfield(Y,'grc_gaba') & strcmp(VM,'all')
                notnan = ~isnan(Y.grc_gaba(:,1,igrc));
                vm = Y.grc_gaba(notnan,:,igrc);
                notdup = ([diff(vm(:,1)); 1] > 0);
                %size(notdup)
                vm = vm(notdup,:);
                Y.GABA(:,igrc) = interp1(vm(:,1),vm(:,2),Y.t);
            end
        end
    end

    % LFP & LFP_G
    if ~isfield(Y,'LFP')
        for k = 1 %:10
            %grc_idx(k).d = find(Y.grc_xyz(:,3)>-50e-6+(k-1)*10e-6 &Y.grc_xyz(:,3)<-50e-6+(k)*10e-6);
            grc_idx(k).d = find(Y.grc_xyz(:,3)>-10e-6 &Y.grc_xyz(:,3)<10e-6)
            %    size(grc_idx(k).d)
        end   % LFP matrix 3Dim: 1 time or vm value, 2 time or vm, 3 cell
        
        for it = 1:length(Y.t) %find(t>50 & t<100)
            if mod(it,100)==0
                disp(it)
            end
            %        for k = 1 %:10
            Y.LFP(:,:,it) = griddata(Y.grc_xyz(grc_idx(1).d,1),Y.grc_xyz(grc_idx(1).d,2),Y.VM(it,grc_idx(1).d),XI,YI);
            if isfield(Y,'GABA') & strcmp(VM,'all')
                Y.LFP_G(:,:,it) = griddata(Y.grc_xyz(grc_idx(1).d,1),Y.grc_xyz(grc_idx(1).d,2),Y.GABA(it,grc_idx(1).d),XI,YI);
            end
        end
    end

    % Make Avi files
    % Open avi file
    af = [Y.dd '/' Y.dd '.avi'];
    mov = avifile(af);

    figure(1)
    for it = 1:size(Y.LFP,3)
        if mod(it,1000)==0
            disp(Y.t(it))
        end
        pcolor(Y.LFP(:,:,it))
        shading interp
        lighting phong        
        %    zlim([vm_min -20])
        caxis([-75 -20])
        colorbar
        title(['Time ' sprintf('%05.3f',Y.t(it)) ' ms'])
        disp(['Time ' sprintf('%05.3f',Y.t(it)) ' ms'])
        %    view(2)
    hold off
    refresh(gcf)
    drawnow update
        pause(0.01)
        %     print('-djpeg90',[Y.dd '/' Y.dd '.jpg'])
        F = getframe(gcf);
        mov = addframe(mov,F);
    end
    mov = close(mov);
    % system(['ffmpeg -y -i ' af ' -target pal-dvd ' af '.avi'])
    % Open avi file
%     af = [Y.dd '/' Y.dd 'GABA.avi'];
%     mov = avifile(af);
% 
%     if isfield(Y,'LFP_G') & strcmp(VM,'all')
%         figure(2)
%         for it = 1:size(Y.LFP_G,3)
%             if mod(it,1000)==0
%                 disp(Y.t(it))
%             end
%             pcolor(Y.LFP_G(:,:,it))
%             shading interp
%             lighting phong
%             %    zlim([vm_min -20])
%             caxis([400 1600])
%             colorbar
%             title(['Time ' sprintf('%05.3f',Y.t(it)) ' ms'])
%             disp(['Time ' sprintf('%05.3f',Y.t(it)) ' ms'])
%             %    view(2)
%             hold off
%             refresh(gcf)
%             drawnow update
%             pause(0.01)
%             %     print('-djpeg90',[Y.dd '/' Y.dd '.jpg'])
%             F = getframe(gcf);
%             mov = addframe(mov,F);
%         end
%         mov = close(mov);
%         % system(['ffmpeg -y -i ' af ' -target pal-dvd ' af '.avi'])
%     end
end