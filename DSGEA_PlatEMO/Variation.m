function [OffDec,OffMask] = Variation(Problem,ParentDec,ParentMask,Fitness,REAL,iter)
% The variation operator for DSGEA (Algorithm 5)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    [N,~]       = size(ParentDec);  
    Parent1Mask = ParentMask(1:floor(end/2),:);
    Parent2Mask = ParentMask(floor(end/2)+1:floor(end/2)*2,:);
    
    [GroupIndex,ChooseGroup,maxGroup] = DynamicGroup(ParentDec,ParentMask,iter);
    
    %% Crossover for mask
    OffMask = Parent1Mask;
    
    chosenindex = GroupIndex == ChooseGroup;
        
    for i = 1 : N/2
        if rand < 0.5            
            index = find(Parent1Mask(i,:)&~Parent2Mask(i,:));
            index = index(TS(-Fitness(index),maxGroup));
            OffMask(i,index) = 0;
        else            
            index = find(~Parent1Mask(i,:)&Parent2Mask(i,:));
            index = index(TS(Fitness(index),maxGroup));
            OffMask(i,index) = Parent2Mask(i,index);
        end
        
    end
   
    
    %% Mutation for mask
    if REAL
        for i = 1 : N/2
            if rand < iter
                index = find(OffMask(i,:)&chosenindex(i,:));
                index = index(TS(-Fitness(index),maxGroup));
                OffMask(i,index) = 0;
            else
                index = find(~OffMask(i,:)&chosenindex(i,:));
                index = index(TS(Fitness(index),maxGroup));
                OffMask(i,index) = 1;
            end
        end      
    end
    
    %% Optimal variation for decision variables
    Parent1Dec  = ParentDec(1:floor(end/2),:);
    Parent2Dec  = ParentDec(floor(end/2)+1:floor(end/2)*2,:);  
    OffDec = DecOptimal(Problem,Parent1Dec,Parent2Dec,GroupIndex,ChooseGroup,REAL);
    
end

function index = TS(Fitness,maxGroup)
% Binary tournament selection helper function
    if isempty(Fitness)
        index = [];
    else
        index = TournamentSelection(maxGroup,1,Fitness);
    end
end