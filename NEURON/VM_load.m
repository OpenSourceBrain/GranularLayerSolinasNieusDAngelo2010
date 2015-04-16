function Y = VM_load(dd,mode)

Y.dd = dd;

Y.grc_xyz = load([dd '/grc_coord.lst']);
Y.goc_xyz = load([dd '/goc_coord.lst']);
Y.glm_xyz = load([dd '/glom_coord.lst']);

Y.goc_targets = textread([dd '/div_goc_targets.lst'],'','emptyvalue',NaN);
Y.glm_target_grcs = textread([dd '/div_glom_target_grc.lst'],'','emptyvalue',NaN);
Y.glm_target_gocs = textread([dd '/div_glom_target_goc.lst'],'','emptyvalue',NaN);
if exist([dd '/Glom_stim.lst'])==2
    Y.glm_stim = load([dd '/Glom_stim.lst']);
    Y.glm_stim = Y.glm_stim+1;
else
    Y.glm_stim = [];
end
% Spike trains
Y.goc_spk = load([dd '/Goc_spks.dat']);
Y.grc_spk = load([dd '/Grc_spks.dat']);
Y.glm_spk = load([dd '/Glom_spks.dat']);

Y.goc = dir([dd '/goc_*_vm.dat']);
Y.grc = dir([dd '/grc_*_vm.dat']);
Y.glm = dir([dd '/glom_*_vm.dat']);

% Find the simulation time length
Y.goc
vm = load([dd '/' Y.goc(1).name]);
stl = max(vm(:,1));
Y.t = (0:0.1:stl)';



if strfind(mode,'V')

    for i = 1:length(Y.goc)
        vm = load([dd '/goc_' mat2str(i-1) '_vm.dat']);
        notdup = ([diff(vm(:,1)); 1] > 0);
        vm = vm(notdup,:);
        % Y.goc_vm(:,i) = interp1(vm(:,1),vm(:,2),Y.t);
        s = size(vm);
        if i>1
            v = size(Y.goc_vm);
            if s(1)>v(1)
                Y.goc_vm = padarray(Y.goc_vm,[s(1)-v(1) 0 0],NaN,'post');
            end
            if s(1)<v(1)
                Y.goc_vm = padarray(Y.goc_vm,[0 0 1],NaN,'post');
            end
            Y.goc_vm(1:s(1),:,i) = vm;
        else
            Y.goc_vm(:,:,i) = vm;
        end
    end
    
    for i = 1:length(Y.grc)
        if mod(i,100)==0
            disp(i)
        end
        vm = load([dd '/grc_' mat2str(i-1) '_vm.dat']);
        notdup = ([diff(vm(:,1)); 1] > 0);
        vm = vm(notdup,:);
        % Y.grc_vm(:,i) = interp1(vm(:,1),vm(:,2),Y.t);
        s = size(vm);
        if i>1
            v = size(Y.grc_vm);
            if s(1)>v(1)
                Y.grc_vm = padarray(Y.grc_vm,[s(1)-v(1) 0 0],NaN,'post');
            end
            if s(1)<v(1)
                Y.grc_vm = padarray(Y.grc_vm,[0 0 1],NaN,'post');
            end
            Y.grc_vm(1:s(1),:,i) = vm;
        else
            Y.grc_vm(:,:,i) = vm;
        end
    end
end

if strfind(mode,'G')
    gb = dir([dd '/Grc_gaba_*_vm.dat']);
    if length(gb)>1
        for i = 1:length(Y.grc)
            if mod(i,100)==0
                disp(i)
            end
            vm = load([dd '/Grc_gaba_' mat2str(i-1) '_vm.dat']);
            notdup = ([diff(vm(:,1)); 1] > 0);
            vm = vm(notdup,:);
            % Y.grc_vm(:,i) = interp1(vm(:,1),vm(:,2),Y.t);
            s = size(vm);
            if i>1
                v = size(Y.grc_gaba);
                if s(1)>v(1)
                    Y.grc_gaba = padarray(Y.grc_gaba,[s(1)-v(1) 0 0],NaN,'post');
                end
                if s(1)<v(1)
                    Y.grc_gaba = padarray(Y.grc_gaba,[0 0 1],NaN,'post');
                end
                Y.grc_gaba(1:s(1),:,i) = vm;
            else
                Y.grc_gaba(:,:,i) = vm;
            end
        end
    end
end
if strfind(mode,'C')
    grc_channels = {'GRANULE_LKG1',...
        'GRANULE_LKG2',...
        'GRANULE_Nmda_leak',...
        'GRANULE_NA',...
        'GRANULE_NAR',...
        'GRANULE_PNA',...
        'GRANULE_KV',...
        'GRANULE_KA',...
        'GRANULE_KIR',...
        'GRANULE_KCA',...
        'GRANULE_KM',...
        'GRANULE_CA'};
    goc_channels = {'Golgi_lkg',...
        'Golgi_Na',...
        'Golgi_NaR',...
        'Golgi_NaP',...
        'Golgi_Ca_HVA',...
        'Golgi_Ca_LVA',...
        'Golgi_KV',...
        'Golgi_KM',...
        'Golgi_KA',...
        'Golgi_BK',...
        'Golgi_SK2',...
        'Golgi_hcn1',...
        'Golgi_hcn2'};
    
    rec_kind = {'_curr','_act','_inact'};
    for k = 1:length(rec_kind), kind = rec_kind(k);
        disp(kind)
        for goc_c = 1:length(goc_channels), chan = goc_channels(goc_c);
            for goc_i = 1:length(Y.goc)
                if mod(goc_i,100)==0
                    disp([kind ' ' mat2str(goc_i)])
                end
                if exist([dd '/goc_' mat2str(goc_i-1) '_' char(chan) char(kind) '.dat'])
                    vm = load([dd '/goc_' mat2str(goc_i-1) '_' char(chan) char(kind) '.dat']);
                    notdup = ([diff(vm(:,1)); 1] > 0);
                    vm = vm(notdup,:);
                    s = size(vm);
                    if goc_i>1
                        v = size(Y.goc_chan(goc_c).goc_kind(k).d);
                        if s(1)>v(1)
                            Y.goc_chan(goc_c).goc_kind(k).d = padarray(Y.goc_chan(goc_c).goc_kind(k).d,[s(1)-v(1) 0 0],NaN,'post');
                        end
                        if s(2)>v(2)
                            Y.goc_chan(goc_c).goc_kind(k).d = padarray(Y.goc_chan(goc_c).goc_kind(k).d,[0 2 0],NaN,'post');
                        end
                        if s(1)<v(1)
                            Y.goc_chan(goc_c).goc_kind(k).d = padarray(Y.goc_chan(goc_c).goc_kind(k).d,[0 0 1],NaN,'post');
                        end
                        size(vm)
                        size(Y.goc_chan(goc_c).goc_kind(k).d)
                        Y.goc_chan(goc_c).goc_kind(k).d(1:s(1),:,goc_i) = vm;
                    else
                        Y.goc_chan(goc_c).goc_kind(k).d(:,:,goc_i) = vm;
                    end
                else
                    Y.goc_chan(goc_c).goc_kind(k).d(1:2,:,goc_i) = nan(2);
                end
            end
        end
    end
    for k = 1:length(rec_kind), kind = rec_kind(k);
        disp(kind)
        for grc_c = 1:length(grc_channels), chan = grc_channels(grc_c);
            for grc_i = 1:length(Y.grc)
                if mod(grc_i,100)==0
                    disp([kind ' ' mat2str(grc_i)])
                end
                if exist([dd '/grc_' mat2str(grc_i-1) '_' char(chan) char(kind) '.dat'])
                    vm = load([dd '/grc_' mat2str(grc_i-1) '_' char(chan) char(kind) '.dat']);
                    notdup = ([diff(vm(:,1)); 1] > 0);
                    vm = vm(notdup,:);
                    s = size(vm);
                    if grc_i>1
                        v = size(Y.grc_chan(grc_c).grc_kind(k).d);
                        if s(1)>v(1)
                            Y.grc_chan(grc_c).grc_kind(k).d = padarray(Y.grc_chan(grc_c).grc_kind(k).d,[s(1)-v(1) 0 0],NaN,'post');
                        end
                        if s(2)>v(2)
                            Y.grc_chan(grc_c).grc_kind(k).d = padarray(Y.grc_chan(grc_c).grc_kind(k).d,[0 2 0],NaN,'post');
                        end
                        if s(1)<v(1)
                            Y.grc_chan(grc_c).grc_kind(k).d = padarray(Y.grc_chan(grc_c).grc_kind(k).d,[0 0 1],NaN,'post');
                        end
                        Y.grc_chan(grc_c).grc_kind(k).d(1:s(1),:,grc_i) = vm;
                    else
                        Y.grc_chan(grc_c).grc_kind(k).d(:,:,grc_i) = vm;
                    end
                else
                    Y.grc_chan(grc_c).grc_kind(k).d(1:2,:,grc_i) = nan(2);
                end
            end
        end
    end
end
