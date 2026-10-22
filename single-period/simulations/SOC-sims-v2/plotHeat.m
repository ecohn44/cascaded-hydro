function plotHeat(x_vals, y_vals, vals, x0, y0, input, xlab, ylab)

close all; 

% Find baseline index
[~, i0] = min(abs(y_vals - y0));
[~, j0] = min(abs(x_vals - x0));

if input == "obj"
    
    J0 = vals(i0,j0);
    Z = 100 * (vals - J0) / abs(J0);
    cblabel = 'Generation Increase (%)';

elseif input == "IVI"
    
    % rescaled to physical parameters
    Z = vals/1e5;
    cblabel = 'Integrated Violation Index (10^5 m^3)';

elseif input == "phi"

    Z = vals;
    cblabel = 'Minimum Reliability Guarantee';

elseif input == "IVI_diff"

    Z = vals;
    cblabel = 'IVI: BON - SSH (m^3)';

elseif input == "energy"

    Z = vals;
    cblabel = 'Energy: SSH - BON (MWh)';
end

figure;
h = imagesc(x_vals, y_vals, Z);
set(h,'AlphaData',~isnan(Z))   % hide NaNs
set(gca,'Color',[0.85 0.85 0.85]) % gray for NaNs

if input == "phi"
    set(gca,'YDir','reverse');
    
    baseline = 0.95;
    d = .3; 
    clim([baseline-d, baseline+d])
    colormap(flipud(parula))

elseif input == "obj"
    set(gca,'YDir','normal');
    colormap(flipud(summer))  

elseif input == "IVI"
    set(gca,'YDir','normal');
    % clim([32, 35.5])
    colormap(winter) 

elseif input == "IVI_diff"
    set(gca,'YDir','normal');
    colormap(parula)  

elseif input == "energy"
    set(gca,'YDir','normal');
    colormap(parula)

else
    set(gca,'YDir','normal');
end

xlabel(xlab,'FontSize',12)
ylabel(ylab,'FontSize',12)

cb = colorbar;
ylabel(cb,cblabel,'FontSize',14)

if input == "phi"
    cb.Limits = [baseline-d baseline];        % zoom legend only
    cb.Ticks  = (baseline-d):0.05:baseline;  
elseif input == "IVI"
    %cb.Limits = [33 35]; 
    %cb.Ticks  = 33:0.5:35;  
end 

grid on
set(gca,'GridColor','k','GridAlpha',0.4,'LineWidth',1.0)

end