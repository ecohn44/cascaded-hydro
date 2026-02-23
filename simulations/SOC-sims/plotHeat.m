function plotHeat(gamma_vals, xi_vals, vals, gamma0, xi0, input)

% Find baseline index
[~, i0] = min(abs(xi_vals - xi0));
[~, j0] = min(abs(gamma_vals - gamma0));

if input == "obj"
    
    % Percent improvement 
    J0 = vals(i0, j0);
    Z = 100 * (vals - J0) / abs(J0);
    cblabel = 'Generation Increase (%)';

elseif input == "IVI"
    
    % Plot IVI as is
    cblabel = 'Integrated Violation Index';
    Z = vals; 
    
end 

figure;
imagesc(gamma_vals, xi_vals, Z);
set(gca,'YDir','normal');

xlabel('Upstream Release Coefficient (\gamma)', 'FontSize', 12);
ylabel('Previous Forecast Error (\xi)', 'FontSize', 12);

cb = colorbar;
ylabel(cb, cblabel, 'FontSize', 14);

grid on;
set(gca,'GridColor','k','GridAlpha',0.4,'LineWidth',1.0);

end