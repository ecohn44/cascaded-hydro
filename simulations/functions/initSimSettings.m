function settings = initSimSettings(season, framework, bounds)
    validSeasons = ["dry", "wet"];
    validFrameworks = ["det", "diu", "ddu"];
    validBounds = ["det", "icc", "jcc-bon", "jcc-ssh"];

    if ~ismember(season, validSeasons)
        error('Invalid season. Choose "DRY" or "WET".');
    end
    if ~ismember(framework, validFrameworks)
        error('Invalid framework. Choose "DET", "DIU", or "DDU".');
    end
    if ~ismember(bounds, validBounds)
        error('Invalid bounds framework".');
    end


    settings = struct('season', season, 'framework', framework, 'bounds', bounds);
end