%%%%%%%%%%%%%%%%%%%%%%% List Of Usefull Variables %%%%%%%%%%%%%%%%%%%%%%%

% empty so far

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% General Opertaions With Data %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%
%%%%%%%%%%%%%%%%%%%%%%% Dataset Paths (Choose One) %%%%%%%%%%%%%%%%%%%%%%%

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A20_d1/experiment_data.h5';
% design 2 !!!

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A20_d2/experiment_data.h5';
% design 1

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A21_d1/experiment_data.h5';
% design 1

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A21_d2/experiment_data.h5';
% design 1

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/1/1_15-15_2_02-04_11-02-47/experiment_data.h5';
% design 1 !!!!!!!!!!!!!!!

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/3/3_15-15_1_02-13_19-41-53/experiment_data.h5';
% design 2

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/3/3_15-15_2_02-14_11-43-05/experiment_data.h5';
% design 1

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/4/4_15-15_1_02-17_12-21-15/experiment_data.h5';
% design 2

%%%%%%%%%%%%%%%%%%%%%%%
% Newest
%%%%%%%%%%%%%%%%%%%%%%%

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A1/a1_p1_d1_03-08_16-51-25/experiment_data.h5';
% design absent

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A1/a1_p2_d1_03-08_17-03-48/experiment_data.h5';
% design 2


% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A2//experiment_data.h5';
% design absent

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A2//experiment_data.h5';
% design absent


% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A3/a3_p1_d1_03-09_17-24-16/experiment_data.h5';
% design 2

path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A3/a3_p1_d2_03-10_19-51-08/experiment_data.h5';
% design 1

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A3/a3_p2_d1_03-09_17-36-05/experiment_data.h5';
% design 1

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A3/a3_p2_d2_03-10_20-02-06/experiment_data.h5';
% design 2


% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A4/a4_p1_d1_03-09_19-02-35/experiment_data.h5';
% design 1

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A4/a4_p1_d2_03-10_18-44-10/experiment_data.h5';
% design 2

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A4/a4_p2_d1_03-09_19-14-09/experiment_data.h5';
% design 2

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/A4/a4_p2_d2_03-10_18-55-09/experiment_data.h5';
% design 1

%%
%%%%%%%%%%%%%%%%%%%%%%% Type Of Design (Choose One) %%%%%%%%%%%%%%%%%%%%%%%

% design 1

Xdl = [2 6 12 16 20 22 26];
Ydl = [4 8 10 14 18 24 28 30];

% design 2

% Xdl = [2 8 12 14 20 22 26];
% Ydl = [4 6 10 16 18 24 28 30];

%%
%%%%%%%%%%%%%%%%%%%%%%% Concatinate Protocoles %%%%%%%%%%%%%%%%%%%%%%%

DMd = hdf5read(path, 'protocol1/raw_data');
sizes(1) = length(DMd);

for i=2:30
    protocolName = sprintf('protocol%d/raw_data', i);
    rom = hdf5read(path, protocolName);
    DMd = [DMd rom];
    sizes(i) = length(rom) + sizes(i-1);
end 

%%
%%%%%%%%%%%%%%%%%%% Protocol Spectrum Visualisation %%%%%%%%%%%%%%%%%%%

% chanal = 23;
% figure

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

%%
%%%%%%%%%%%%%%%%%%%%%%% Total Spectrum Visualisation %%%%%%%%%%%%%%%%%%%%%%%

% chanal = 23;
% figure
% 
% Xc = hdf5read(path, 'protocol2/raw_data');
% Yc = hdf5read(path, 'protocol4/raw_data');
% 
% for i=1:7
%     protocolNameXc = sprintf('protocol%d/raw_data', Xdl(i));
%     protocolNameYc = sprintf('protocol%d/raw_data', Ydl(i));
%     Xc = [Xc hdf5read(path, protocolNameXc)];
%     Yc = [Yc hdf5read(path, protocolNameYc)];
% end
% 
% [PxxX,WX] = pwelch(Xc(chanal, :), 3000);
% [PxxY,WY] = pwelch(Yc(chanal, :), 3000);
% subplot(1, 2, 1)
% plot(WX/pi*500/2, PxxX)
% title('Total Real Power')
% xlim([0 50])
% ylim([0 10*10^-10])
% subplot(1, 2, 2)
% plot(WY/pi*500/2, PxxY)
% title('Total Mock Power')
% xlim([0 50])
% ylim([0 10*10^-10])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Permutation Test (F(x)=CSP(Real, Mock)) for different frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% General Settings %%%%%%%%%%%%%%%%%%%%%%%

steps = 300; % number of permutations
srate = 500; % sampling rate of data
h = waitbar(0, 'Wait...'); % initializing the process bar
protTotNum = 4; % number of protocols to be considered
N_sub_prot = 4; % Number of subprotocols. A protocol length 
                % will be divided by this number. 
distDel = 500; % number of possible distorted by a filter values that 
               % should be droped
example = hdf5read(path, 'protocol2/raw_data');
protocol_lenght = length(example);
tic; % start timer

%%%%%%%%%%%%%%%%%%%%%%% Permutation Test %%%%%%%%%%%%%%%%%%%%%%%

for frequency=1:40 % frequency for FIR filter
    
    % Setting FIR filter 
    band=[frequency frequency+2];
    order = 400;
    [b]= fir1(order, band*2/srate);
    a = [1];

    % Appling FIR filter to dataset
    
    DM = filtfilt(b, a, DMd')';
    
    % Choose one of the following options. Comment the other ones. 
    
    % 1. Protocol Analysis
    
%     N_prot = 4;
%     
%     for i=1:N_prot
%        protLength = 7504;
%        sliseX = sizes(Xdl(i)-1):(sizes(Xdl(i)-1)+protLength);
%        sliseY = sizes(Ydl(i)-1):(sizes(Ydl(i)-1)+protLength);
%        X{i} = DM(:, sliseX);
%        Y{i} = DM(:, sliseY);
%     end
    
    % 2. Sub-Protocol Analysis
    
%         N_sub_prot = 4;
%         protLength = 7504;
%         
%         protocolX = 1;
%         vectorX = sizes(Xdl(protocolX)-1):(sizes(Xdl(protocolX)-1)+protLength);
%         sliseX = DM(:, vectorX);
%         lenX = round(length(sliseX)/N_sub_prot);
%         for i=1:N_sub_prot
%             indexX = (lenX*(i-1)+1):(lenX*i);
%             X{i} = sliseX(:, indexX);
%         end
% 
%         protocolY = 1;
%         vectorY = sizes(Ydl(protocolY)-1):(sizes(Ydl(protocolY)-1)+protLength);
%         sliseY = DM(:, vectorY);
%         lenY = round(length(sliseY)/N_sub_prot);
%         for i=1:N_sub_prot
%             indexY = (lenY*(i-1)+1):lenY*i;
%             Y{i} = sliseY(:, indexY);
%         end

    % 2. Mixed Analysis
    
    counter = 0;
    
    for prot=1:protTotNum
            protLength = protocol_lenght; % Length of the Real-Feedback or 
                                % Mock-Feedback protocol, that will be used
                                % in the analysis.

            protocolX = prot; % Number of the Real-Feedback protocol,
                              % contained in Xdl vector.
            protocolY = prot; % Number of the Mock-Feedback protocol,
                              % contained in Ydl vector.
                              
            vectorX = sizes(Xdl(protocolX)-1):(sizes(Xdl(protocolX)-1)+protLength);
            % Number of samples, that will be taken to make a slise from
            % common data matrix to get a particular protocol. 
            vectorY = sizes(Ydl(protocolY)-1):(sizes(Ydl(protocolY)-1)+protLength);
            % Number of samples, that will be taken to make a slise from
            % common data matrix to get a particular protocol. 
            
            sliseX = DM(:, vectorX); % Get the slise with values got in the
                                     % previous step. 
            sliseY = DM(:, vectorY); % Get the slise with values got in the
                                     % previous step.
            
            lenX = floor(length(sliseX)/N_sub_prot);
            lenY = floor(length(sliseY)/N_sub_prot);
            
            for i=1:N_sub_prot
                indexX = (lenX*(i-1)+1):(lenX*i);
                indexY = (lenY*(i-1)+1):(lenY*i);
                X{i+counter} = sliseX(:, indexX);
                Y{i+counter} = sliseY(:, indexY);
            end
            
            counter = counter + N_sub_prot;
    end

    % Making permutations
    
    for step=1:steps

        ind = randperm(2*counter);
        for cell_maker = 1:counter
            XY{cell_maker} = X{cell_maker};
            XY{cell_maker+counter} = Y{cell_maker};
        end
    
        XX = [XY{ind(1)}]; YY = [XY{ind(counter+1)}];
        
        for vector_maker = 2:counter
            XX = [XX XY{ind(vector_maker)}];
            YY = [YY XY{ind(counter+vector_maker)}];
        end

        l_XX = length(XX);
        l_YY = length(YY);

        % find covariances C10 and C20 
        C10 = XX*XX'/l_XX;
        C20 = YY*YY'/l_YY;

        nchan = size(C10,1);

        C10n = C10/trace(C10);
        C20n = C20/trace(C20);

        % Tikhonov regularization
        C1 = C10n + 0.05 * trace(C10n) * eye(nchan) / size(C10n,1); 
        C2 = C20n + 0.05 * trace(C20n) * eye(nchan) / size(C20n,1);
        % try different regularization parameters p_reg

        [V, d] = eig(C1,C2);
        D_sg(step,:) = diag(d)'; 
        
        % Getting a 3D matrix of components (x) for different frequencies
        % for different permutations (z)
        
        for ii=1:size(d,2)
           D_perm(frequency,ii,step) = d(ii,ii);
        end
        
        waitbar((((frequency-1)*300)+step)/(40*300)) 
    end

    % Getting real eigenvalues for different frequencies
    
    ind = 1:counter*2;
    
    for cell_maker=1:counter
        XYr{cell_maker} = X{cell_maker};
        XYr{cell_maker+counter} = Y{cell_maker};
    end
    
    Xr = [XYr{ind(1)}]; Yr = [XYr{ind(counter+1)}];
    
    for vector_maker = 2:counter
        Xr = [Xr XYr{ind(vector_maker)}];
        Yr = [Yr XYr{ind(counter+vector_maker)}];
    end
    
    l_Xr = length(Xr);
    l_Yr = length(Yr);

    % Find covariances C10 and C20 
    C10r = Xr*Xr'/l_Xr;
    C20r = Yr*Yr'/l_Yr;

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
    Dr(frequency, :) = diag(de)';
    Wr{frequency}= inv(Ve');
    FF(frequency) = frequency+1;
    
    for numCompStat=1:30
        
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making Normalization And New P-value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make Vector Of Maximum

% for step = 1:steps
%     DmxDistr(step) = max(max(D(:,:,step)));
% end;

% Make New P-Value

% for f= 1:40
%     for numCompStat=1:30
%         comp = De(f, numCompStat);
%         p_value_norm(f, numCompStat) = length(find(DmxDistr > comp)) / steps;
%     end
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chanlocs_vis = makeChanlocsVis();

H = uicontrol('Style', 'PushButton', ...
                    'String', 'Break', ...
                    'Callback', 'delete(gcbf)');         
    
% Plot "component's eigenvalue = f(frequency)"

figure(1);
imagesc(p_valuePlus);
title('Right Sided');colorbar;

% figure(2);
% imagesc(p_valuePlus);
% title('Right Sided');
% colorbar;
% figure(3);
% imagesc(p_valueMinus);
% title('Left Sided');
% colorbar;


while (ishandle(H))
    
    figure(1);
    
    % Picking figure's coordinates
    [x0,y0] = ginput(1);
    x = round(x0);
    y = round(y0);
    
    if x > length(Wr)
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nameDir = '/Users/basilminkov/Desktop/Neurofeedback/Analysis/a3_p1_d2_03-10_19-51-08(half)';

if exist(nameDir) ~= 7
   mkdir(nameDir)
end
save('/Users/basilminkov/Desktop/Neurofeedback/Analysis/a3_p1_d2_03-10_19-51-08(half)/a3_p1_d2_03-10_19-51-08.mat')

% saveas(figure(1), '/Users/basilminkov/Desktop/Neurofeedback/Analysis/A21_d2/A21_d2.jpg');
