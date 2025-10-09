function settings = initSimSettings(season, method, framework, bounds)
    validSeasons = ["dry", "wet"];
    validMethods = ["minlp", "pwl"];
    validFrameworks = ["det", "diu", "ddu"];
    validBounds = ["det", "icc", "jcc-bon", "jcc-ssh"];

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

    settings = struct('season', season, 'method', method, 'framework', framework, 'bounds', bounds);
end