%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting permutated values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X{1} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A12_d1', 'protocol2'); 
% X{2} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A12_d1', 'protocol8'); 
% X{3} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A12_d1', 'protocol12');
% X{4} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A12_d1', 'protocol14');
% Y{1} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A12_d1', 'protocol4'); 
% Y{2} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A12_d1', 'protocol6'); 
% Y{3} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A12_d1', 'protocol10'); 
% Y{4} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A12_d1', 'protocol16');

X{1} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A18_d1\raw.h5', 'protocol2'); 
X{2} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A18_d1\raw.h5', 'protocol8'); 
X{3} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A18_d1\raw.h5', 'protocol12');
X{4} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A18_d1\raw.h5', 'protocol14');
Y{1} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A18_d1\raw.h5', 'protocol4'); 
Y{2} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A18_d1\raw.h5', 'protocol6'); 
Y{3} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A18_d1\raw.h5', 'protocol10'); 
Y{4} = hdf5read('D:\Minkov\Alpha Feedback\Data\One Day Data\A18_d1\raw.h5', 'protocol16');

h = waitbar(0, 'Wait...');
steps = 100;
srate = 500;

for f = 1:40

    band=[f f+2];
    order = 3;
    [b]= fir1(90, band*2/srate);
    a = [1];
    
    for step=1:steps

        ind = randperm(8);
        XY = {X{1} X{2} X{3} X{4} Y{1} Y{2} Y{3} Y{4}};

        XXd = [XY{ind(1)} XY{ind(2)} XY{ind(3)} XY{ind(4)}];
        YYd = [XY{ind(5)} XY{ind(6)} XY{ind(7)} XY{ind(8)}];
        XXd(18,:) = []; % subject's special
        YYd(18,:) = []; % subject's special

        XX = filtfilt(b,a,XXd')';
        YY = filtfilt(b,a,YYd')';

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
        W{f, step} = inv(V');
        % time_passed = toc

        D_sg(step,:) = diag(d)';
        
        for ii=1:size(d,2)
           D(f,ii,step) = d(ii,ii);
        end
    end %step
   
    D_f{f} = D_sg
    
    ind = 1:8;
    XY = {X{1} X{2} X{3} X{4} Y{1} Y{2} Y{3} Y{4}};

    Xd = [XY{ind(1)} XY{ind(2)} XY{ind(3)} XY{ind(4)}];
    Yd = [XY{ind(5)} XY{ind(6)} XY{ind(7)} XY{ind(8)}];
    Xd(18,:) = [];
    Yd(18,:) = [];
    XX = Xd;
    YY = Yd;

    l_X = length(XX);
    l_Y = length(YY);

    % Find covariances C10 and C20 
    C10 = XX*XX'/l_X;
    C20 = YY*YY'/l_Y;

    nchan = size(C10,1);

    C10n = C10/trace(C10);
    C20n = C20/trace(C20);

    % Tikhonov regularization
    C1 = C10n + 0.1 * trace(C10n) * eye(nchan) / size(C10n,1); 
    C2 = C20n + 0.1 * trace(C20n) * eye(nchan) / size(C20n,1);
    % try different regularization parameters p_reg

    % W0 - cells of matrixes of eigenvectors
    % D0 - a matrix of eigenvalues: component (rows), frequency (column)
    % FF - vector of frequencies
    [Ve, de] = eig(C1,C2); 
    De(f, :) = diag(de)';
    W0{f}= inv(Ve');
    FF(f) = f+1;
    
    for numCompStat=1:30
        
        for i=1:steps
            permEVs(i) = D_sg(i, numCompStat);
        end

        comp = De(f, numCompStat);
        
        p_value_lol(f, numCompStat) = 1 - ((1*comp)/max(permEVs));
        p_value(f, numCompStat) = (length(find(permEVs>comp))/length(permEVs));
        
    end
    waitbar((f/40)) 
end

close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making Normalization And New P-value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make Vector Of Maximum

for step = 1:steps
    DmxDistr(step) = max(max(D(:,:,step)));
end;

% Make New P-Value

for f= 1:40
    for numCompStat=1:30
        comp = De(f, numCompStat);
        p_value_norm(f, numCompStat) = length(find(DmxDistr > comp)) / steps;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Some Stuff 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = p_value;

% Setting topoplotting properties
chanlocs_vis = makeChanlocsVis();

H = uicontrol('Style', 'PushButton', ...
                    'String', 'Break', ...
                    'Callback', 'delete(gcbf)');         
    
% Plot "component's eigenvalue = f(frequency)"

f1 = figure(1);
% plot(FF, D0(:, numComp));
imagesc(p);
colorbar;

while (ishandle(H))
    
    f1;
    
    % Picking figure's coordinates
    [x0,y0] = ginput(1);
    x = round(x0);
    y = round(y0);
    
    if x > length(W0)
        continue
    end
    
    if x < 0
        continue
    end
    
    if y > length(W0)
        continue
    end
    
    if y < 0
        continue
    end
    
    fprintf('Fq: %d', y)
    fprintf('Comp: %d', x)
    
    figure(2);
    clf();
    topoplot(W0{y}(:, x), chanlocs_vis, 'electrodes', 'on');
    title(num2str(De(y, x)));
    
end

saveas(f1, 'D:\Minkov\Alpha Feedback\Analysis\A18_d1\A18_d1_p_value.png')
save('D:\Minkov\Alpha Feedback\Analysis\A18_d1\A18_d1.mat')
