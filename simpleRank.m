function [finalRank,scores] = simpleRank( data )
%SIMPLERANKBYCLASS Summary of this function goes here
%   Detailed explanation goes here

% Data format: vector of structures containing the fields 'data' and
% 'label', being that label is an integer in [1, numLabel] and data is a 
% (row) vector of real numbers

    nEx = length(data);

    labels = zeros(nEx,1);
    for thisEx = 1 : nEx
        labels(thisEx) = data(thisEx).label;
    end
    nLabel = length(unique(labels));

    display('Performing NN to calculate scores');
    % Find scores
    [scores] = nn_LOO_SimpleRank(data,nLabel);
    
    display('Ranking the examples for each class label');
    rankForLabel = [];
    for thisLabel = 1 : nLabel
        
        exOfThisLabel = find(labels == thisLabel);
        [~,thisRank] = sort(scores(exOfThisLabel), 'descend');
        
        rankForLabel = [rankForLabel; struct('pos', 1, ...
            'rank', thisRank, 'classIdx', exOfThisLabel)];
    end
    

    display('Constructing the final rank');
    finalRank = zeros(nEx,1);
    pos = 1;

    while pos <= nEx
        
        % Get one example of each (available) class
        thisRound = [];
        for thisLabel = 1 : nLabel

            thisStr = rankForLabel(thisLabel);
            if (thisStr.pos <= length(thisStr.rank))
                
                thisRound = [thisRound; ...
                    thisStr.classIdx(thisStr.rank(thisStr.pos))];
                
                rankForLabel(thisLabel).pos = ...
                    rankForLabel(thisLabel).pos + 1;
                
            end

        end
        
        % sort the examples by the score
        [~,thisRoundRank] = sort(scores(thisRound), 'descend');
        
        finalRank(pos:pos+length(thisRound)-1) = ...
            thisRound(thisRoundRank);
        pos = pos+length(thisRound);

                
    end

end

