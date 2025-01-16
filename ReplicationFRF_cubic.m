set(findobj(0,'type','figure'),'visible','on')
close all

clearvars;
clc;


%<<<<<<<<<<<<<<<<<<<<<<<<<< Part 1. Import Data >>>>>>>>>>>>>>>>>>>>>>>>>>%
frf_dataset_import = importdata('pb_d_joint_table.csv');
frf_data = frf_dataset_import.data;
dN = (length(frf_dataset_import.colheaders)-1)/2;   % dimension of cross-sections

for j = 1:length(frf_dataset_import.colheaders)
    eval([frf_dataset_import.colheaders{j} '= frf_data(:,j);'])
end

start_date = 1;
end_date = 58;    % 73 if the final year is 2022; 58 if the final year is 2007

balance  = frf_data(start_date+1:end_date,2:dN+1);          % y variable: primary balance-to-GDP ratio
debt = frf_data(start_date:end_date-1,end-dN+1:end);        % x variable: public debt-to-GDP ratio

countries = cell(dN,1);              
countries{1} = 'Austria';    
countries{2} = 'Germany';              
countries{3} = 'Norway'; 
countries{4} = 'Portugal'; 
countries{5} = 'Switzerland';   

% countries{3} = 'Netherlands';
% countries{4} = 'Norway'; 
% countries{5} = 'Portugal'; 
% countries{6} = 'Switzerland'; 


% Fully modified estimation + joint Bonferroni KPSS tests + Wald statistics
Vhat = diff(debt,1);
vx = debt(2:end,:)';
dT = size(Vhat,1);               % sample size due to one lag

ind_selected = 1:dN;
mZ = zeros(dN*dT,4*dN);
for t = 1:dT
    mZ(dN*(t-1)+1:dN*t,:) = repmat([ones(dN,1) vx(ind_selected',t) vx(ind_selected',t).^2 vx(ind_selected', t).^3],1,dN).*kron(eye(dN),ones(1,4));
end
mY = balance(:,ind_selected)'; 
mY = mY(:,2:end);
vy = mY(:);

[vbetaSOLS,vbetaSUR,vbetaGLS,KpssJoint,...
    WaldstatLinear,WaldstatQuadra,WaldstatCubic,...
    Table] = fm_inference_cubic(dN,dT,vy,mZ,vx,Vhat,...
    'OrganizeOutputs','YES','SelectedCountries',countries(ind_selected)');

% Check significance level: the number indicates how many stars, 3 stars
% for alpha = 1%
sig_level = [0.1; 0.05; 0.01];
pvalue_linear = zeros(dN,3);
pvalue_quadra = zeros(dN,3);
pvalue_cubic = zeros(dN,3);
for id = 1:length(sig_level)
    dAlpha = sig_level(id);
    pvalue_linear(WaldstatLinear > chi2inv(1-dAlpha,1)) = id;
    pvalue_quadra(WaldstatQuadra > chi2inv(1-dAlpha,1)) = id;
    pvalue_cubic(WaldstatCubic> chi2inv(1-dAlpha,1)) = id;
end
if end_date == 73
    writetable(Table,'FullEstCubic22.xls');
elseif end_date == 58
    writetable(Table,'FullEstCubic07.xls');
else
    disp('Check the starting and ending years')
end




% Plots
% vuhat_SOLS = vy-mZ*vbetaSOLS;
% mUhat_SOLS = reshape(vuhat_SOLS,[dN,dT]);

mYfit_SOLS = reshape(mZ*vbetaSOLS,[dN,dT]);
mYfit_SUR = reshape(mZ*vbetaSUR,[dN,dT]);
mYfit_GLS = reshape(mZ*vbetaGLS,[dN,dT]);
for id = 1:dN
    figure (id)
    plot(year(start_date+2:end_date), mY(id,:), ':k','LineWidth',3)
    hold on 
    plot(year(start_date+2:end_date), mYfit_SOLS(id,:), '--b','LineWidth',4)
    plot(year(start_date+2:end_date), mYfit_SUR(id,:), '-.m','LineWidth',4)
    plot(year(start_date+2:end_date), mYfit_GLS(id,:), '-*r','LineWidth',4)
    hold off
    grid minor
    lg = legend({'$pb$','$\widehat{pb}_{\mathrm{SOLS}}$',...
        '$\widehat{pb}_{\mathrm{SUR}}$','$\widehat{pb}_{\mathrm{GLS}}$'},...
        'Interpreter','latex','Location','Best',...
        'FontSize',35,'Orientation','horizontal');
    set(lg,'color','none','Box','off')
    xlim([year(start_date+2),year(end_date)])
    % ylim([-3,4])
    xlabel('year','Interpreter','latex')
    if end_date == 73
        plot_title = '1950 -- 2022';
        figure_name = '_cubic22';
    elseif end_date == 58
        plot_title = '1950 -- 2007';
        figure_name = '_cubic07';
    else
        disp('Check the starting and ending years')
    end
    title(plot_title,'FontSize',40,'Interpreter','latex')
    ax = gca;
    ax.FontSize = 30;
    xaxisproperties= get(ax, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex';
    yaxisproperties= get(ax, 'YAxis');
    yaxisproperties.TickLabelInterpreter = 'latex';
    set(gcf,'Position',[0,0,600,600])
    saveas(gcf,[countries{id},figure_name,'.png'])
end





