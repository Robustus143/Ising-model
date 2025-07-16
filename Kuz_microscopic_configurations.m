
% Generate microscopic spin configurations at different times

L = 100;
N = L^2;
sweepCount = 20000;
saveSteps = 2000:2000:18000; % 9 snapshots
temps = [2.0, 2.5];

% Color map: -1 -> black, +1 -> blue
customCMap = [0 0 0; 0.3 0.6 1];

for Tidx = 1:2
    T = temps(Tidx);
    beta = 1/T;
    config = sign(rand(L)-0.5); % random init
    snapList = {};
    stepIndex = 1;

    for sweep = 1:sweepCount
        config = metropolisStep(config, beta);
        if ismember(sweep, saveSteps)
            snapList{stepIndex} = config;
            stepIndex = stepIndex + 1;
        end
    end

    figure('Name', sprintf('Microscopic Configurations at T = %.1f', T), 'Color', 'w');
    for i = 1:9
        subplot(3,3,i);
        imagesc(snapList{i});
        axis square off;
        colormap(customCMap);
        caxis([-1 1]);
        title(sprintf('t = %d', saveSteps(i)));
    end

    sgtitle(sprintf('Microscopic configurations for T = %.1f', T));
end

% --- Metropolis step function ---
function config = metropolisStep(config, beta)
    L = size(config,1);
    for k = 1:L^2
        i = randi(L); j = randi(L);
        s = config(i,j);
        nb = config(mod(i,L)+1,j) + config(mod(i-2,L)+1,j) + ...
             config(i,mod(j,L)+1) + config(i,mod(j-2,L)+1);
        dE = 2 * s * nb;
        if dE < 0 || rand < exp(-beta * dE)
            config(i,j) = -s;
        end
    end
end
