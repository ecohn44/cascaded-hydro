function settings = initSimSettings(season, framework, bounds, reference)
    validSeasons = ["dry", "wet"];
    validFrameworks = ["det", "diu", "ddu"];
    validBounds = ["det", "icc", "jcc-bon", "jcc-ssh"];
    validRef = ["mean", "envelope"];

    if ~ismember(season, validSeasons)
        error('Invalid season. Choose "DRY" or "WET".');
    end
    if ~ismember(framework, validFrameworks)
        error('Invalid framework. Choose "DET", "DIU", or "DDU".');
    end
    if ~ismember(bounds, validBounds)
        error('Invalid bounds framework".');
    end
    if ~ismember(reference, validRef)
        error('Invalid tracking reference framework".');
    end

    settings = struct('season', season, 'framework', framework, 'bounds', bounds, 'ref', reference);
end