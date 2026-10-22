clear; clc; close all;

load("resultsBonferroni/monteCarloResults.mat","results");

%% Case to inspect
selected_year  = 2018;
selected_theta = 5;
selected_kappa = 1;

years  = 2018:2024;
thetas = [0 5 10 15];

y = find(years == selected_year,1);
h = find(thetas == selected_theta,1);
k = find(results.kappa == selected_kappa,1);

m_diu = find(results.frameworks == "diu",1);
m_ddu = find(results.frameworks == "ddu",1);

T       = size(results.V,2);
n_units = size(results.V,1);
t       = 1:T;

% Average every trajectory across Monte Carlo scenarios (dimension 7)
V_diu   = reshape(mean(results.V(:,:,y,h,m_diu,k,:),7),n_units,T);
V_ddu   = reshape(mean(results.V(:,:,y,h,m_ddu,k,:),7),n_units,T);
p_diu   = reshape(mean(results.p(:,:,y,h,m_diu,k,:),7),n_units,T);
p_ddu   = reshape(mean(results.p(:,:,y,h,m_ddu,k,:),7),n_units,T);
u_diu   = reshape(mean(results.u(:,:,y,h,m_diu,k,:),7),n_units,T);
u_ddu   = reshape(mean(results.u(:,:,y,h,m_ddu,k,:),7),n_units,T);
sp_diu  = reshape(mean(results.sp(:,:,y,h,m_diu,k,:),7),n_units,T);
sp_ddu  = reshape(mean(results.sp(:,:,y,h,m_ddu,k,:),7),n_units,T);
q_diu   = reshape(mean(results.q(:,:,y,h,m_diu,k,:),7),n_units,T);
q_ddu   = reshape(mean(results.q(:,:,y,h,m_ddu,k,:),7),n_units,T);
std_diu = reshape(mean(results.std(:,:,y,h,m_diu,k,:),7),n_units,T);
std_ddu = reshape(mean(results.std(:,:,y,h,m_ddu,k,:),7),n_units,T);

SOC_mean = results.SOC_mean(:,:,y);
SOC_p10  = results.SOC_p10(:,:,y);
SOC_p90  = results.SOC_p90(:,:,y);

diu_color = [0.00 0.45 0.74];
ddu_color = [0.85 0.33 0.10];

for i = 1:n_units
    sys = results.sysparams(i);

    % Forebay elevation calculated from the mean volume trajectory
    Vnorm_diu = (V_diu(i,:) - sys.min_V)/(sys.max_V - sys.min_V);
    Vnorm_ddu = (V_ddu(i,:) - sys.min_V)/(sys.max_V - sys.min_V);
    Vnorm_diu = min(1,max(0,Vnorm_diu));
    Vnorm_ddu = min(1,max(0,Vnorm_ddu));
    head_diu = sys.min_h + (sys.max_h - sys.min_h).*Vnorm_diu.^sys.b;
    head_ddu = sys.min_h + (sys.max_h - sys.min_h).*Vnorm_ddu.^sys.b;

    figure("Color","w","Position",[100 100 1200 600]);
    tl = tiledlayout(2,3,"TileSpacing","compact","Padding","compact");

    nexttile;
    plot(t,u_diu(i,:),"Color",diu_color,"LineWidth",1.5); hold on;
    plot(t,u_ddu(i,:),"Color",ddu_color,"LineWidth",1.5);
    yline(sys.min_ut,"r--","HandleVisibility","off");
    yline(sys.max_ut,"r--","HandleVisibility","off");
    title("Generation Outflow"); ylabel("Flow"); grid on;
    legend("DIU","DDU","Location","best");

    nexttile;
    plot(t,p_diu(i,:),"Color",diu_color,"LineWidth",1.5); hold on;
    plot(t,p_ddu(i,:),"Color",ddu_color,"LineWidth",1.5);
    yline(sys.F,"r--","HandleVisibility","off");
    title("Hydropower Generation"); ylabel("Power"); grid on;
    legend("DIU","DDU","Location","best");

    nexttile;
    plot(t,sp_diu(i,:),"Color",diu_color,"LineWidth",1.5); hold on;
    plot(t,sp_ddu(i,:),"Color",ddu_color,"LineWidth",1.5);
    title("Spill Outflow"); ylabel("Flow"); grid on;
    legend("DIU","DDU","Location","best");

    nexttile;
    plot(t,V_diu(i,:),"Color",diu_color,"LineWidth",1.5); hold on;
    plot(t,V_ddu(i,:),"Color",ddu_color,"LineWidth",1.5);
    plot(t,SOC_mean(i,:),"k--","LineWidth",1.2);
    plot(t,SOC_p10(i,:),"Color",[0.60 0.60 0.60],"LineStyle",":");
    plot(t,SOC_p90(i,:),"Color",[0.60 0.60 0.60],"LineStyle",":");
    yline(sys.min_V,"r--","HandleVisibility","off");
    yline(sys.max_V,"r--","HandleVisibility","off");
    title("Reservoir Volume"); ylabel("Volume"); grid on;
    legend("DIU","DDU","SOC mean","SOC p10","SOC p90", ...
        "Location","best");
    
    nexttile;
    yyaxis left
    plot(t,q_diu(i,:),"Color",diu_color,"LineWidth",1.5); hold on;
    plot(t,q_ddu(i,:),"Color","black","LineWidth",1.5);
    ylabel("Realized Inflow");

    yyaxis right
    plot(t,std_diu(i,:),"--","Color",diu_color,"LineWidth",1.2);
    plot(t,std_ddu(i,:),"Color",ddu_color,"LineWidth",1.2);
    ylabel("Forecast Standard Deviation");
    title("Inflow and Forecast Uncertainty"); grid on;
    legend("DIU inflow","DDU inflow","DIU std","DDU std", ...
        "Location","best");

    nexttile;
    plot(t,head_diu,"Color",diu_color,"LineWidth",1.5); hold on;
    plot(t,head_ddu,"Color",ddu_color,"LineWidth",1.5);
    yline(sys.min_h,"r--","HandleVisibility","off");
    yline(sys.max_h,"r--","HandleVisibility","off");
    title("Forebay Elevation"); ylabel("Elevation"); grid on;
    legend("DIU","DDU","Location","best");

    axs = findall(gcf,"Type","axes");
    set(axs,"XLim",[1 T],"FontName","Times New Roman","FontSize",9);
    for ax = axs'
        xlabel(ax,"Hour");
    end

    title(tl,sprintf("Unit %d: year %d, theta = %g, kappa = %g", ...
        sys.unit,selected_year,selected_theta,selected_kappa), ...
        "FontWeight","bold");
end
