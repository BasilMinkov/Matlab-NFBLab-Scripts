%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting real & permutated values for different frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reading all the protocols from the dataset. I'm planing to make a few
% functions for each experimental design in order not to rewrite this part
% for each particepent.
% Xn and Yn ??protocols for different conditions. 
% "Xd" ? "dirty" X dataset. "Yd" ? "dirty" Y dataset. 

Xd{1} = hdf5read('/Users/basilminkov/Desktop/Neurofeedback/Data/A12/raw.h5', 'protocol2'); 
Xd{2} = hdf5read('/Users/basilminkov/Desktop/Neurofeedback/Data/A12/raw.h5', 'protocol8'); 
Xd{3} = hdf5read('/Users/basilminkov/Desktop/Neurofeedback/Data/A12/raw.h5', 'protocol12');
Xd{4} = hdf5read('/Users/basilminkov/Desktop/Neurofeedback/Data/A12/raw.h5', 'protocol14');
Yd{1} = hdf5read('/Users/basilminkov/Desktop/Neurofeedback/Data/A12/raw.h5', 'protocol4'); 
Yd{2} = hdf5read('/Users/basilminkov/Desktop/Neurofeedback/Data/A12/raw.h5', 'protocol6'); 
Yd{3} = hdf5read('/Users/basilminkov/Desktop/Neurofeedback/Data/A12/raw.h5', 'protocol10'); 
Yd{4} = hdf5read('/Users/basilminkov/Desktop/Neurofeedback/Data/A12/raw.h5', 'protocol16');

% General Settings

steps = 300; % number of permutations
srate = 500; % sampling rate of data
h = waitbar(0, 'Wait...'); % initializing the process bar
distDel = 500; % number of possible distorted values that should be droped


for frequency=1:40 % frequency for FIR filter

    % Setting FIR filter 
    band=[frequency frequency+2];
    order = 600;
    [b]= fir1(order, band*2/srate);
    a = [1];

    % Appling FIR filter to two kinds of datasets 
    
    for i = 1:length(Xd)
        x = filtfilt(b, a, Xd{i}')';
        x(18,:) = []; % subject's special bad channel
        X{i} = x(:, distDel:end-distDel); % getting rid of distortion
        y = filtfilt(b, a, Yd{i}')';
        y(18,:) = []; % subject's special bad channel
        Y{i} = y(:, distDel:end-distDel); % getting rid of distortion
    end 
    
    % Making permutations
    
    for step=1:steps

        ind = randperm(8);
        XY = {X{1} X{2} X{3} X{4} Y{1} Y{2} Y{3} Y{4}};

        XX = [XY{ind(1)} XY{ind(2)} XY{ind(3)} XY{ind(4)}];
        YY = [XY{ind(5)} XY{ind(6)} XY{ind(7)} XY{ind(8)}];

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
    
    ind = 1:8;
    XYr = {X{1} X{2} X{3} X{4} Y{1} Y{2} Y{3} Y{4}};

    Xr = [XYr{ind(1)} XYr{ind(2)} XYr{ind(3)} XYr{ind(4)}];
    Yr = [XYr{ind(5)} XYr{ind(6)} XYr{ind(7)} XYr{ind(8)}];
    
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

p = p_valueMinus;

% Setting topoplotting properties
chanlocs_vis = makeChanlocsVis();

% H = uicontrol('Style', 'PushButton', ...
%                     'String', 'Break', ...
%                     'Callback', 'delete(gcbf)');         
    
% Plot "component's eigenvalue = f(frequency)"

figure(1);
% plot(FF, D0(:, numComp));
imagesc(p);
colorbar;

% while (ishandle(H))
%     
%     figure(1);
%     
%     % Picking figure's coordinates
%     [x0,y0] = ginput(1);
%     x = round(x0);
%     y = round(y0);
%     
%     if x > length(W0)
%         continue
%     end
%     
%     if x < 0
%         continue
%     end
%     
%     if y > length(W0)
%         continue
%     end
%     
%     if y < 0
%         continue
%     end
%     
%     fprintf('Fq: %d', y)
%     fprintf('Comp: %d', x)
%     
%     figure(2);
%     clf();
%     topoplot(W0{y}(:, x), chanlocs_vis, 'electrodes', 'on');
%     title(num2str(De(y, x)));
%     
% end

% save('/Users/basilminkov/Desktop/Neurofeedback/Analysis/A12/A12_explicit/A12_explicit_p_value_lol', p_value_lol);
% save('/Users/basilminkov/Desktop/Neurofeedback/Analysis/A12/A12_explicit/A12_explicit_p_value', p_value);
% save('/Users/basilminkov/Desktop/Neurofeedback/Analysis/A12/A12_explicit/A12_explicit_W0', W0);
% save('/Users/basilminkov/Desktop/Neurofeedback/Analysis/A12/A12_explicit/A12_explicit_De', De);
% saveas('/Users/basilminkov/Desktop/Neurofeedback/Analysis/A12/A12_explicit/A12_explicit_comp=f(Fr)', figure(1));
