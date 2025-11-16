function [Fitness] = FitnessUpdate(Fitness,Population,iter)
% Fitness update strategy for DSGEA (Algorithm 3)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Update fitness based on population sparsity
    t = ceil(mean(sum(logical(Population.decs),2)));
    
    Fitness = max(1,Fitness-t*iter*mode(logical(Population.decs)));
    
    Fitness = Fitness+t*iter*mode(~logical(Population.decs));
end