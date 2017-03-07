%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting real & permutated values for different frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/3/3_15-15_1_02-13_19-41-53/experiment_data.h5';
% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/3/3_15-15_2_02-14_11-43-05/experiment_data.h5';
% path = '/Users/basilminkov/Desktop/Neurofeedback/Data/1/1_15-15_2_02-04_11-02-47/experiment_data.h5';
path = '/Users/basilminkov/Desktop/Neurofeedback/Data/4/4_15-15_1_02-17_12-21-15/experiment_data.h5';

% General Settings

steps = 300; % number of permutations
srate = 500; % sampling rate of data
h = waitbar(0, 'Wait...'); % initializing the process bar
distDel = 500; % number of possible distorted values that should be droped

% design 1

% Xdl = [2 6 12 16 20 22 26];
% Ydl = [4 8 10 14 18 24 28 30];

% design 2

Xdl = [2 8 12 14 20 22 26];
Ydl = [4 6 10 16 18 24 28 30];

DMd = hdf5read(path, 'protocol1/raw_data');
sizes(1) = length(DMd);

for i=2:30
    protocolName = sprintf('protocol%d/raw_data', i);
    rom = hdf5read(path, protocolName);
    DMd = [DMd rom];
    sizes(i) = length(rom) + sizes(i-1);
end 

for frequency=1:40 % frequency for FIR filter
    
    % Setting FIR filter 
    band=[frequency frequency+2];
    order = 400;
    [b]= fir1(order, band*2/srate);
    a = [1];

    % Appling FIR filter to two kinds of datasets 
    
    DM = filtfilt(b, a, DMd')';
    
    % Choose one of the following options. Comment the other one. 
    
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
    counter = 0
    
    for prot=1:7
            N_sub_prot = 4;
            protLength = 7500;

            protocolX = prot;
            vectorX = sizes(Xdl(protocolX)-1):(sizes(Xdl(protocolX)-1)+protLength);
            sliseX = DM(:, vectorX);
            lenX = round(length(sliseX)/N_sub_prot);
            for i=1:N_sub_prot
                indexX = (lenX*(i-1)+1):(lenX*i);
                X{i+counter} = sliseX(:, indexX);
            end

            protocolY = prot;
            vectorY = sizes(Ydl(protocolY)-1):(sizes(Ydl(protocolY)-1)+protLength);
            sliseY = DM(:, vectorY);
            lenY = round(length(sliseY)/N_sub_prot);
            for i=1:N_sub_prot
                indexY = (lenY*(i-1)+1):lenY*i;
                Y{i+counter} = sliseY(:, indexY);
            end
            
            counter = counter + N_sub_prot;
    end

    % Making permutations
    
    for step=1:steps

        ind = randperm(2*counter);
        for cell_maker = 1:counter
            XY{cell_maker} = X{cell_maker};
            XY{cell_maker*2} = Y{cell_maker};
        end
    
        XX = [XY{ind(1)}]; YY = [XY{ind(16+1)}];
        
        for vector_maker = 2:counter
            XX = [XX XY{ind(vector_maker)}];
            YY = [YY XY{ind(16+vector_maker)}];
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
        
    end

    % Getting real eigenvalues for different frequencies
    
    ind = 1:counter*2;
    
    for cell_maker=1:counter
        XYr{cell_maker} = X{cell_maker};
        XYr{cell_maker*2} = Y{cell_maker};
    end
    
    Xr = [XYr{ind(1)}]; Yr = [XYr{ind(16+1)}];
    
    for vector_maker = 2:counter
        Xr = [Xr XYr{ind(vector_maker)}];
        Yr = [Yr XYr{ind(16+vector_maker)}];
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
        p_valueDouble(frequency, numCompStat) = (length(find(permEVs<comp))/length(permEVs));

    end
    
    waitbar((frequency/40)) 
    
end

close(h);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Some Stuff 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 
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


% mkdir('/Users/basilminkov/Desktop/Neurofeedback/Analysis/1_15-15_2')
% % saveas(figure(1), '/Users/basilminkov/Desktop/Neurofeedback/Analysis/A21_d2/A21_d2.jpg');
% save('/Users/basilminkov/Desktop/Neurofeedback/Analysis/1_15-15_2/1_15-15_2.mat')

