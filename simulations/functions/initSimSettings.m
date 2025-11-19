function settings = initSimSettings(season, drought, method, framework, bounds)
    validSeasons = ["dry", "wet"];
    validMethods = ["minlp", "pwl"];
    validFrameworks = ["det", "diu", "ddu"];
    validBounds = ["det", "icc", "jcc-bon", "jcc-ssh"];
    validDroughts = ["constant", "pulse", "extended"];

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

    settings = struct('season', season, 'drought', drought, 'method', method, 'framework', framework, 'bounds', bounds);
end