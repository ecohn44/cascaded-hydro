function [a_nb, d_nb] = nearestNeighbor(phi_bon, a, d)
    [rows, cols] = find(~isnan(phi_bon));
    if isempty(rows)
        error('No converged BON cells found.');
    end
    distances = sqrt((rows - a).^2 + (cols - d).^2);
    [~, idx]  = min(distances);
    a_nb = rows(idx);
    d_nb = cols(idx);
end