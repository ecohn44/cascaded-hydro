function settings = initSimSettings(season, drought, method, framework, bounds, volPrice)
    validSeasons = ["dry", "wet"];
    validMethods = ["minlp", "pwl"];
    validFrameworks = ["det", "diu", "ddu"];
    validBounds = ["det", "icc", "jcc-bon", "jcc-ssh"];
    validDroughts = ["constant", "pulse", "extended"];
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
    if ~ismember(drought, validDroughts)
        error('Invalid drought framework".');
    end
    if ~ismember(volPrice, validVolumePrices)
        error('Invalid volume pricing framework".');
    end

    settings = struct('season', season, 'drought', drought, 'method', method, 'framework', framework, 'bounds', bounds, 'volPrice', volPrice);
end