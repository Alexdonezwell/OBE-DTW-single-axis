function plot_adding_random(dsName, ratioVec, seeds)
%PLOT_ADDING_RANDOM Summary of this function goes here
%   Detailed explanation goes here

    nRuns = length(seeds);
    ratioVec = sort(ratioVec);
    
    acc_rates_dtw = zeros(length(ratioVec),1);
    acc_rates_cdtw = zeros(length(ratioVec),1);
    acc_rates_obedtw = zeros(length(ratioVec),1);
    acc_rates_obe_cdtw = zeros(length(ratioVec),1);
    acc_rates_psidtw = zeros(length(ratioVec),1);
    acc_rates_cdtw_psidtw = zeros(length(ratioVec),1);
    
    for i = 1 : length(ratioVec)
        
        for j = 1 : nRuns
        
            load(['exp_seed', num2str(seeds(j)), '/', ...
                dsName, '_', num2str(ratioVec(i)), ...
                '_results.mat'], 'acc_dtw','acc_cdtw','acc_obedtw','acc_obe_cdtw',  ...
                'acc_psidtw', 'acc_psicdtw');

            acc_rates_dtw(i) = acc_rates_dtw(i) + acc_dtw;
            
            acc_rates_cdtw(i) = acc_rates_cdtw(i) + acc_cdtw;
            
            acc_rates_obedtw(i) = acc_rates_obedtw(i) + acc_obedtw;
            
            acc_rates_obe_cdtw(i) = acc_rates_obe_cdtw(i) + acc_obe_cdtw;
            
            acc_rates_psidtw(i) = acc_rates_psidtw(i) + acc_psidtw;
            
            acc_rates_cdtw_psidtw(i) = acc_rates_cdtw_psidtw(i) ...
                + acc_psicdtw;
            
        end
        
    end
    
    % average
    acc_rates_dtw = acc_rates_dtw / nRuns;
    acc_rates_cdtw = acc_rates_cdtw / nRuns;
    acc_rates_obedtw = acc_rates_obedtw / nRuns;
    acc_rates_obe_cdtw = acc_rates_obe_cdtw / nRuns;
    acc_rates_psidtw = acc_rates_psidtw / nRuns;
    acc_rates_cdtw_psidtw = acc_rates_cdtw_psidtw / nRuns;

    
    FigHandle = handle(figure);
    set(FigHandle, 'Position', [250 250 700 300]);
  
    plot(acc_rates_dtw, 'k', 'LineWidth', 3.25); 
    hold on;
    plot(acc_rates_cdtw, 'LineWidth', 2.25, 'Color', [0 0.75 0]);
    plot(acc_rates_obedtw, 'm', 'LineWidth', 1.25);
    plot(acc_rates_obe_cdtw, 'c', 'LineWidth', 2.75);
    plot(acc_rates_cdtw_psidtw, 'r', 'LineWidth', 1.5);
    plot(acc_rates_psidtw, 'b', 'LineWidth', 1);

    xlabel('Relative length of added random walk');
    ylabel('Accuracy');
    
    NumTicks = length(ratioVec);
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    set(gca,'XTickLabel',ratioVec)
    set(gca,'box','off')
    legend('Location','Best', 'DTW','cDTW','OBEDTW','OBE-cDTW',  ...
        'psi-cDTW', 'psi-DTW');
    
end

