function settings = initSimSettings(season, method, framework)
    validSeasons = ["dry", "wet"];
    validMethods = ["minlp", "pwl"];
    validFrameworks = ["det", "diu", "ddu"];

    if ~ismember(season, validSeasons)
        error('Invalid season. Choose "DRY" or "WET".');
    end
    if ~ismember(method, validMethods)
        error('Invalid method. Choose "MINLP" or "PWL".');
    end
    if ~ismember(framework, validFrameworks)
        error('Invalid framework. Choose "DET", "DIU", or "DDU".');
    end

    settings = struct('season', season, 'method', method, 'framework', framework);
end