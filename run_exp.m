function run_exp(dsName, ratioVec, seed)

    for i = 1 : length(ratioVec)

        load(['exp_seed', num2str(seed), '/', dsName, '_', ...
            num2str(ratioVec(i)),'.mat'])
        
        display('OBEDTW')
        [acc_obedtw] = nnOBEDTW(finalData,0.1);
        
        display('OBE-cDTW')
        [acc_obe_cdtw] = nnOBE_cDTW(finalData,0.2,0.1);

        display('PSI-cDTW')
        [acc_psicdtw] = nnPSIcDTW(finalData, 0.1);

        display('PSI-DTW')
        [acc_psidtw] = nnPSIDTW(finalData, 0.1);

        display('cDTW')
        [acc_cdtw] = nncDTW(finalData, 0.1);
        
        display('DTW')
        [acc_dtw] = nnDTW(finalData);


        
        save(['exp_seed', num2str(seed), '/', dsName, '_', ...
            num2str(ratioVec(i)), '_results.mat'], 'acc_psicdtw', ...
            'acc_psidtw','acc_cdtw',  'acc_dtw', 'acc_obedtw','acc_obe_cdtw');

    end
    
    clear; clc;
  
end

