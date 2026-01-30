function settings = initSimSettings(season, scenario, method, framework, bounds, volPrice)
    validSeasons = ["dry", "wet"];
    validMethods = ["minlp", "pwl"];
    validFrameworks = ["det", "diu", "ddu"];
    validBounds = ["det", "icc", "jcc-bon", "jcc-ssh"];
    validScenarios = ["constant", "pulse", "extended"];
    validVolumePrices = ["none", "static", "dynamic"];

    if ~ismember(season, validSeasons)
        error('Invalid season. Choose "DRY" or "WET".');
    end
    if ~ismember(method, validMethods)
        error('Invalid method. Choose "MINLP" or "PWL".');
    end
    if ~ismember(framework, validFrameworks)
        error('Invalid framework. Choose "DET", "DIU", or "DDU".');
    end
    if ~ismember(bounds, validBounds)
        error('Invalid bounds framework".');
    end
    if ~ismember(scenario, validScenarios)
        error('Invalid drought framework".');
    end
    if ~ismember(volPrice, validVolumePrices)
        error('Invalid volume pricing framework".');
    end

    settings = struct('season', season, 'scenario', scenario, 'method', method, 'framework', framework, 'bounds', bounds, 'volPrice', volPrice);
end