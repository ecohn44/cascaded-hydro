clear; clc; close all;

load("resultsBonferroni/monteCarloResults.mat","results");

frameworks = results.frameworks;
kappa = results.kappa;
T = size(results.V,2);
M = length(frameworks);

%% Number of hours completed by each scenario
hours_run = T*ones(size(results.failed));
hours_run(results.failed) = results.failure_time(results.failed);

%% 1. Physical-bound failure rate
failure_rate = 100*mean(results.failed,3);

%% 2. Historical-envelope IVI per hour
IVI_rate = results.IVI./hours_run;
mean_IVI = mean(IVI_rate,3);

%% 3. Generation per hour
generation = squeeze(sum(results.p,[1 2],"omitnan"));
generation_rate = generation./hours_run;
mean_generation = mean(generation_rate,3,"omitnan");

%% Plot
figure;
tiledlayout(1,3);

metrics = {failure_rate,mean_IVI,mean_generation};
titles = {"Physical Failure Rate","Envelope IVI per Hour", ...
          "Generation per Hour"};
ylabels = {"Failure Rate (%)","Normalized IVI","Generation"};

for j = 1:3
    nexttile;
    hold on;

    for m = 1:M
        plot(kappa,metrics{j}(m,:),"-o","LineWidth",1.5);
    end

    xlabel("\kappa");
    ylabel(ylabels{j});
    title(titles{j});
    grid on;
end

legend(frameworks,"Location","eastoutside");