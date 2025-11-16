function [GroupIndex,ChooseGroup,maxGroup] = DynamicGroup(ParentDec,ParentMask,iter)
% Dynamic grouping strategy for DSGEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate max group size
    maxGroup = max(2,ceil(iter*(mean(sum(ParentMask,2))*0.5*size(ParentMask,1)/size(ParentMask,2))));
    
    %% Perform dynamic grouping
    GroupIndex = [];
    noOfSolutions = 0.5*size(ParentMask,1);
    numberOfVariables = size(ParentMask,2);
    
    for sol = 1:noOfSolutions
        varsPerGroup = floor(numberOfVariables/maxGroup);
        vars = sum(ParentMask);
        [~,I] = sort(vars);
        
        outIndexList = ones(1,numberOfVariables);
        for i = 1:maxGroup-1
            % Grouping based on sorted variable indices
            outIndexList(I(((i-1)*varsPerGroup)+1:i*varsPerGroup)) = i;
        end
        outIndexList(I(((maxGroup-1)*varsPerGroup)+1:end)) = maxGroup;
        GroupIndex = [GroupIndex;outIndexList];
    end
    
    ChooseGroup = randi(maxGroup,size(GroupIndex,1),1);
end