% 2D Ising Model Simulation using Metropolis Monte Carlo Algorithm
% Author: Volodymyr Kuz
% Course: Kinetic Theory and Stochastic Simulations

%% Simulation Parameters
latticeSize = 100;
totalSpins = latticeSize^2;
sweepCount = 20000;
temperatureSet = [2.0, 2.5];
temperatureGrid = linspace(0.1, 4.0, 40);

rng(1); % Fix seed for reproducibility

%% Utility Functions

function spinMatrix = initSpinConfiguration(sizeL, mode)
    if mode == 1
        spinMatrix = ones(sizeL);
    elseif mode == 2
        spinMatrix = -ones(sizeL);
    else
        spinMatrix = sign(rand(sizeL) - 0.5);
    end
end

function spinMatrix = applyMetropolis(spinMatrix, invTemp)
    L = size(spinMatrix,1);
    for step = 1:L^2
        x = randi(L); y = randi(L);
        currentSpin = spinMatrix(x,y);
        neighbours = spinMatrix(mod(x,L)+1,y) + spinMatrix(mod(x-2,L)+1,y) + ...
                     spinMatrix(x,mod(y,L)+1) + spinMatrix(x,mod(y-2,L)+1);
        deltaE = 2 * currentSpin * neighbours;
        if deltaE < 0 || rand < exp(-invTemp * deltaE)
            spinMatrix(x,y) = -currentSpin;
        end
    end
end

function eVal = calculateEnergy(spinMatrix)
    L = size(spinMatrix,1); eVal = 0;
    for row = 1:L
        for col = 1:L
            s = spinMatrix(row,col);
            nbs = spinMatrix(mod(row,L)+1,col) + spinMatrix(row,mod(col,L)+1);
            eVal = eVal - s * nbs;
        end
    end
end

function mVal = calculateMagnetization(spinMatrix)
    mVal = sum(spinMatrix(:));
end

%% Time Evolution at Fixed Temperatures

for T = temperatureSet
    beta = 1/T;
    figure('Name',sprintf('Temporal Evolution at T = %.1f',T));
    tiledlayout(1,2);

    for startCase = 1:3
        repeatNum = 2 * (startCase == 3) + 1 * (startCase < 3);
        for instance = 1:repeatNum
            config = initSpinConfiguration(latticeSize, startCase);
            magnetHistory = zeros(1, sweepCount);
            energyHistory = zeros(1, sweepCount);
            for t = 1:sweepCount
                config = applyMetropolis(config, beta);
                magnetHistory(t) = calculateMagnetization(config)/totalSpins;
                energyHistory(t) = calculateEnergy(config)/totalSpins;
            end
            nexttile(1)
            plot(magnetHistory, 'DisplayName', sprintf('Init%d-%d',startCase,instance)); hold on
            title(sprintf('Magnetization at T = %.1f', T)); xlabel('Steps'); ylabel('Magnet./Spin')
            nexttile(2)
            plot(energyHistory, 'DisplayName', sprintf('Init%d-%d',startCase,instance)); hold on
            title(sprintf('Energy at T = %.1f', T)); xlabel('Steps'); ylabel('Energy/Spin')
        end
    end
    nexttile(1); legend; nexttile(2); legend;
end

%% Thermal Averages Across Temperature Grid
avgMag = zeros(size(temperatureGrid));
avgEnergy = zeros(size(temperatureGrid));
chiVals = zeros(size(temperatureGrid));
heatVals = zeros(size(temperatureGrid));

for idx = 1:length(temperatureGrid)
    T = temperatureGrid(idx); beta = 1/T;
    mSamples = []; eSamples = [];
    for trial = 1:8
        config = initSpinConfiguration(latticeSize, 3);
        for s = 1:sweepCount
            config = applyMetropolis(config, beta);
            if s > sweepCount/2
                m = calculateMagnetization(config)/totalSpins;
                e = calculateEnergy(config)/totalSpins;
                mSamples(end+1) = m;
                eSamples(end+1) = e;
            end
        end
    end
    avgMag(idx) = movmean(abs(mean(mSamples)), 5);
    avgEnergy(idx) = mean(eSamples);
    chiVals(idx) = beta*totalSpins*(mean(mSamples.^2) - mean(mSamples)^2);
    heatVals(idx) = beta^2*totalSpins*(mean(eSamples.^2) - mean(eSamples)^2);
end

% exact result for mean magnetization below Tc
Tc = 2 / log(1 + sqrt(2)); % ≈ 2.269
exactM = zeros(size(temperatureGrid));
for i = 1:length(temperatureGrid)
    T = temperatureGrid(i);
    if T < Tc
        exactM(i) = (1 - sinh(2/T)^(-4))^(1/8);
    else
        exactM(i) = 0;
    end
end

figure('Name','Analytical vs Simulated Results');
subplot(2,2,1);
plot(temperatureGrid, avgMag, 'bo-', 'DisplayName', 'MC Simulation'); hold on;
plot(temperatureGrid, exactM, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Exact Result');
xlabel('Temperature'); ylabel('<m>'); title('Magnetization');
legend;

subplot(2,2,2);
plot(temperatureGrid, avgEnergy, 'bo-');
xlabel('Temperature'); ylabel('<e>'); title('Energy');

subplot(2,2,3);
plot(temperatureGrid, chiVals, 'bo-');
xlabel('Temperature'); ylabel('\chi'); title('Susceptibility');

subplot(2,2,4);
plot(temperatureGrid, heatVals, 'bo-');
xlabel('Temperature'); ylabel('C'); title('Specific Heat');
