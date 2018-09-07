function create_ds_adding_random( dsDir, dsName, ratioVec, seed )
% Loads the dataset (referred by dsName)
% Adds (controlled) random walk after and before the series
% The rw length is half the size of the larger ratio referred by ratioVec
% For each value in ratioVec, cuts the series and save a z-norm dataset

    if (nargin < 4)
        seed = 1;
    end
    rng(seed);
    
    if (~exist(['exp_seed', num2str(seed)], 'dir'))
        mkdir(['exp_seed', num2str(seed)]);
    end
        
    if (dsDir(end) ~= '\' && dsDir(end) ~= '/')
        dsDir = [dsDir, '\'];
    end

    if (max(ratioVec) > 1 || min(ratioVec) < 0)
        display('All values in ratioVec (3rd param) shall be in [0,1]');
        return;
    end
     
              trainData = load(['/Users/apple/Desktop/UCR_code/UCR_TS_Archive_2015/Cricket_X/Cricket_X_TRAIN' ]);
              testData = load(['/Users/apple/Desktop/UCR_code/UCR_TS_Archive_2015/Cricket_X/Cricket_X_TEST']);
     
%           trainData = load(['/Users/apple/Desktop/UCR_code_prune/UCR_TS_Archive_2015/uWaveGestureLibrary_X/uWaveGestureLibrary_X_TRAIN' ]);
%           testData = load(['/Users/apple/Desktop/UCR_code_prune/UCR_TS_Archive_2015/uWaveGestureLibrary_X/uWaveGestureLibrary_X_TEST']);
%               trainData = csvread('/Users/apple/Desktop/UCR_code/UCR_TS_Archive_2015/Cricket_X/feeding_z_training.csv');
%               testData = csvread('/Users/apple/Desktop/UCR_code/UCR_TS_Archive_2015/Cricket_X/feeding_z_testing.csv');

    
    trLabels = trainData(:,1);
    trainData = trainData(:,2:end);
    
    teLabels = testData(:,1);
    testData = testData(:,2:end);
    
    origLen = size(trainData,2);
        
    maxRWLen = ceil(max(ratioVec * origLen));
    
    % Adding RW to the training data     
    trainData = addRW(trainData, maxRWLen, seed);
    
    % Adding RW to the test data     
    testData = addRW(testData, maxRWLen);
    
    for i = 1 : length(ratioVec)
        
        thisLen = origLen + ceil(ratioVec(i) * origLen);
        midPoint = size(testData,2)/2;

        iniPoint = ceil(midPoint - thisLen/2 + 1);
        endPoint = floor(midPoint + thisLen/2);
        
        % Data to save
        trData = znorm(trainData(:,iniPoint:endPoint));
        teData = znorm(testData(:,iniPoint:endPoint));
        
        finalTeData = [];
        for j = 1 : length(teLabels)
            finalTeData = [finalTeData; struct('data', teData(j,:), ...
                'label', teLabels(j))];
        end
        
        finalTrData = [];
        for j = 1 : length(trLabels)
            finalTrData = [finalTrData; struct('data', trData(j,:), ...
                'label', trLabels(j))];
        end
        finalData = struct( 'train', finalTrData, 'test', finalTeData);
        
        save(['exp_seed', num2str(seed), '/', ...
            dsName, '_', num2str(ratioVec(i)), '.mat'], ...
            'finalData');
        
    end

end




function [ts] = znorm(ts)

    if (size(ts,1) > 1)
        for i = 1 : size(ts,1)
            ts(i,:) = znorm(ts(i,:));
        end
    else
        ts = (ts - mean(ts)) ./ std(ts);
    
    end
end