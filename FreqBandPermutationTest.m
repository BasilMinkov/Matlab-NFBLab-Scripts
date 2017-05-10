%%%%%%%%%%%%%%%%%%%%%%% List Of Usefull Variables %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% General Opertaions With Data %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

warning('Check channels, N_sub_prot, data_object.channels_list!')

%%%%%%%%%%%%%%%%%%%%%%% XML Object %%%%%%%%%%%%%%%%%%%%%%%

data_object = eegData();
data_object.path = '/Users/basilminkov/Desktop/Neurofeedback/data_30min/a42_d1_03-24_18-52-01/';
data_object = data_object.makeParsing();
data_object.channels_list(1:end-2)
% [usefulProtocolsList, numberProtocolList, numberList, encodedProtocolList] = data_object.getUsefulProtocolsList();

%%%%%%%%%%%%%%%%%%%%%%% Concatinate Protocoles %%%%%%%%%%%%%%%%%%%%%%%

% should be inside of the class

DMd = [];
id = 0;
indices(1) = 1;
for i=1:length(data_object.protocols_list)
    ram_protocol = hdf5read([data_object.path data_object.h5_filename], ['protocol' int2str(i) '/raw_data']);
    id = id + length(ram_protocol); 
    indices(i+1) = id;
    DMd = [DMd ram_protocol];
end

clear i id ram_protocol 

%%%%%%%%%%%%%%%%%%% Protocol Spectrum Visualisation %%%%%%%%%%%%%%%%%%%

% chanal = 23;
% figure
% 
% for i=1:7
%     protocolNameXc = sprintf('protocol%d/raw_data', Xdl(i));
%     protocolNameYc = sprintf('protocol%d/raw_data', Ydl(i));
%     Xc = hdf5read(path, protocolNameXc);
%     Yc = hdf5read(path, protocolNameYc);
%     [PxxX,WX] = pwelch(Xc(chanal, :), 3000);
%     [PxxY,WY] = pwelch(Yc(chanal, :), 3000);
%     subplot(7, 2, (-1)+2*i)
%     plot(WX/pi*500/2, PxxX)
%     title(sprintf('Protocol %d Real', Xdl(i)))
%     xlim([0 50])
%     ylim([0 10*10^-10])
%     subplot(7, 2, 2*i)
%     plot(WY/pi*500/2, PxxY)
%     title(sprintf('Protocol %d Mock', Ydl(i)))
%     xlim([0 50])
%     ylim([0 10*10^-10])
% 
% 
% end

%%%%%%%%%%%%%%%%%%%%%%% ? Spectrum Visualisation %%%%%%%%%%%%%%%%%%%%%%%

% chanal = 15;
% figure
% 
% Xc = hdf5read([data_object.path data_object.h5_filename], 'protocol2/raw_data');
% Yc = hdf5read([data_object.path data_object.h5_filename], 'protocol3/raw_data');
% 
% [PxxX,WX] = pwelch(Xc(chanal, :), 3000);
% [PxxY,WY] = pwelch(Yc(chanal, :), 3000);
% subplot(1, 2, 1)
% plot(WX/pi*500/2, PxxX)
% title('? Real Power')
% xlim([0 50])
% % ylim([0 10*10^-10])
% subplot(1, 2, 2)
% plot(WY/pi*500/2, PxxY)
% title('? Mock Power')
% xlim([0 50])
% % ylim([0 10*10^-10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Permutation Test (F(x)=CSP(Real, Mock)) for different frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% General Settings %%%%%%%%%%%%%%%%%%%%%%%

min_frequency = 1; 
frequencies = 30;
channels = 32;
steps = 300; % number of permutations
srate = 500; % sampling rate of data
h = waitbar(0, 'Wait...'); % initializing the process bar
N_sub_prot = 3; % Number of subprotocols. A protocol length 
                % will be divided by this number. 
distDel = 500; % number of possible distorted by a filter values that 
               % should be droped
tic; % start timer

%%%%%%%%%%%%%%%%%%%%%%% Prepare Indices %%%%%%%%%%%%%%%%%%%%%%%

% should be inside of the class

for i = 1:length(data_object.protocols_list)
    if regexp(data_object.protocols_list{i}, 'FBR') > 0 % 'FB0''FBR'
        RNs(i) = i;
    end                    
end

for i = 1:length(data_object.protocols_list)
    if regexp(data_object.protocols_list{i}, 'FBM\d+') > 0 % 'FBM\d+' 'FB[1-9]'
        MNs(i) = i;
    end                    
end

RNs = RNs(RNs~=0);
MNs = MNs(MNs~=0);

Rids = indices(RNs);
Mids = indices(MNs);

counterR = 1;
for i = 1:length(Rids)
    protocol_length = indices(RNs(i)+1) - indices(RNs(i));
    subprot_lenght = floor(protocol_length/N_sub_prot);
    for j = 1:N_sub_prot
        new_Rids(counterR) = Rids(i) + (j-1)*subprot_lenght;
        counterR = counterR + 1;
    end
end

counterM = 1;
for i = 1:length(Mids)
    protocol_length = indices(MNs(i)+1) - indices(MNs(i));
    subprot_lenght = floor(protocol_length/N_sub_prot);
    for j = 1:N_sub_prot
        new_Mids(counterM) = Mids(i) + (j-1)*subprot_lenght;
        counterM = counterM + 1;
    end
end

Rids = new_Rids;
Mids = new_Mids;
counterR = counterR - 1;
counterM = counterM - 1;

% Rids = Rids(counterR/2+20:end);
% Mids = Mids(counterM/2+20:end);
% counterR = length(Rids);
% counterM = length(Mids);

counterG = counterR+counterM;
Gids = [Rids Mids];

D_sg = zeros(steps, channels);
Dr = zeros(frequencies, channels);
FF = zeros(frequencies, 1);

%%%%%%%%%%%%%%%%%%%%%%% Permutation Test %%%%%%%%%%%%%%%%%%%%%%%

for frequency=min_frequency:frequencies % frequency for FIR filter
    
    % Setting FIR filter 
    band=[frequency frequency+2];
    order = 400;
    b = fir1(order, band*2/srate);
    a = 1;

    % Appling FIR filter to dataset
    
    DM = filtfilt(b, a, DMd')';
    
    % Making permutations
    
    for step=1:steps
        
%         tic
        % so far subprot_lenght are equal for all the protocols
        
        inds = randperm(counterG);
        
        ind = inds(1:counterG/2);
        R = zeros(channels, length(ind)*subprot_lenght);
        for i = 1:length(ind)
            R(:, 1+(i-1)*subprot_lenght:i*subprot_lenght) = DM(:, Gids(ind(i)):Gids(ind(i))+(subprot_lenght-1));
        end
        
        ind = inds((counterG/2)+1:counterG);
        M = zeros(channels, length(ind)*subprot_lenght);
        for i = 1:length(ind)
            M(:, 1+(i-1)*subprot_lenght:i*subprot_lenght) = DM(:, Gids(ind(i)):Gids(ind(i))+(subprot_lenght-1));
        end
                
        l_R = length(R);
        l_M = length(M);

        % find covariances C10 and C20 
        C10 = R*R'/l_R;
        C20 = M*M'/l_M;

        nchan = size(C10,1);

        C10n = C10/trace(C10);
        C20n = C20/trace(C20);

        % Tikhonov regularization
        C1 = C10n + 0.05 * trace(C10n) * eye(nchan) / size(C10n,1); 
        C2 = C20n + 0.05 * trace(C20n) * eye(nchan) / size(C20n,1);
        % try different regularization parameters p_reg

        [V, d] = eig(C1,C2); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%         d(1, 1)
%         toc
        D_sg(step,:) = diag(d)'; 
        
        waitbar((((frequency-1)*steps)+step)/(frequencies*steps)) 
    end
    
    % Getting real eigenvalues
    
    ind = 1:counterR;
    Rr = zeros(channels, counterR*subprot_lenght);
    for i = 1:counterR
        Rr(:, 1+(i-1)*subprot_lenght:i*subprot_lenght) = DM(:, Rids(ind(i)):Rids(ind(i))+(subprot_lenght-1));
    end

    ind = 1:counterM;
    Mr = zeros(channels, counterM*subprot_lenght);
    for i = 1:counterM
        Mr(:, 1+(i-1)*subprot_lenght:i*subprot_lenght) = DM(:, Mids(ind(i)):Mids(ind(i))+(subprot_lenght-1));
    end
    
    l_Rr = length(Rr);
    l_Mr = length(Mr);

    % Find covariances C10 and C20 
    C10r = Rr*Rr'/l_Rr;
    C20r = Mr*Mr'/l_Mr;

    nchan = size(C10r,1);

    C10nr = C10r/trace(C10r);
    C20nr = C20r/trace(C20r);

    % Tikhonov regularization
    C1r = C10nr + 0.1 * trace(C10nr) * eye(nchan) / size(C10nr,1); 
    C2r = C20nr + 0.1 * trace(C20nr) * eye(nchan) / size(C20nr,1);
    % try different regularization parameters p_reg

    % Wr - cells of matrixes of eigenvectors
    % Dr - a matrix of eigenvalues: component (rows), frequency (column)
    % FF - vector of frequencies
    % ...r means "real"
    [Ve, de] = eig(C1r,C2r); 
%     de(1,1)
    
    Dr(frequency, :) = diag(de)';
    Wr{frequency} = inv(Ve');
    FF(frequency) = frequency+1;
    
    for numCompStat=1:channels

        for i=1:steps
            permEVs(i) = D_sg(i, numCompStat);
        end

        comp = Dr(frequency, numCompStat);
        
        p_valuePlus(frequency, numCompStat) = (length(find(permEVs>comp))/length(permEVs));
        p_valueMinus(frequency, numCompStat) = (length(find(permEVs<comp))/length(permEVs));

    end
    
    disp(sprintf('frequency %d is done!', frequency))
    toc   
end

close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ind = 1:counterR;
% Rr = zeros(32, counterR*subprot_lenght);
% for i = 1:counterR
%     Rr(:, 1+(i-1)*subprot_lenght:i*subprot_lenght) = DM(:, Rids(ind(i)):Rids(ind(i))+(subprot_lenght-1));
% end
% 
% ind = 1:counterM;
% Mr = zeros(32, counterM*subprot_lenght);
% for i = 1:counterM
%     Mr(:, 1+(i-1)*subprot_lenght:i*subprot_lenght) = DM(:, Mids(ind(i)):Mids(ind(i))+(subprot_lenght-1));
% end
% 
% [icaR, W, T, mu] = fastICA(Rr,32);
% [icaM, W, T, mu] = fastICA(Mr,32);
% 
% for i = 1:32
%     [PxxX,WX] = pwelch(icaR(i, :), 3000);
%     subplot(4, 8, i)
%     plot(WX/pi*500/2, PxxX)
%     title(['Real component ' int2str(i)])
%     xlim([0 50])
%     % ylim([0 10*10^-10])
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nameDir = data_object.path;

if exist(nameDir) ~= 7
   mkdir(nameDir)
end
save([data_object.path 'matlab_workspace.mat'], 'p_valuePlus', 'Dr', 'Wr', 'min_frequency')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data_object.channels_list = data_object.channels_list(1:end-2)
chanlocs_vis = makeChanlocsVis(data_object.channels_list);

H = uicontrol('Style', 'PushButton', ...
                    'String', 'Break', ...
                    'Callback', 'delete(gcbf)');         
    
% Plot "component's eigenvalue = f(frequency)"

figure(1);
imagesc(p_valuePlus);
title('Right Sided');colorbar;

while (ishandle(H))
    
    figure(1);
    
    % Picking figure's coordinates
    [x0,y0] = ginput(1);
    x = round(x0);
    y = round(y0);
    
    if x > length(Wr{min_frequency})
        continue
    end
    
    if x < 0
        continue
    end
    
    if y > length(Wr)
        continue
    end
    
    if y < 0
        continue
    end
    
    fprintf('Fq: %d', y)
    fprintf('Comp: %d', x)
    
    figure(2);
    clf();
    topoplot(Wr{y}(:, x), chanlocs_vis, 'electrodes', 'on');
    title(num2str(Dr(y, x)));
    
end
