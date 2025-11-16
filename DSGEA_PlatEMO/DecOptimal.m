function [Offspring] = DecOptimal(Problem,Parent1,Parent2,GroupIndex,ChooseGroup,REAL)
% Optimal variation for decision variables in DSGEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %% Parameter setting
    [proC,disC,~,disM] = deal(1,20,1,20);
    [N,D]   = size(Parent1);

    if REAL~=1
        %% Genetic operators for binary encoding
        Offspring = ones(N,D);
    else
        %% Genetic operators for real encoding (SBX and Polynomial Mutation)
        
        % Simulated Binary Crossover (SBX)
        beta = zeros(N,D);
        mu   = rand(N,D);
        beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
        beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
        beta = beta.*(-1).^randi([0,1],N,D);
        beta(rand(N,D)<0.5) = 1;
        beta(repmat(rand(N,1)>proC,1,D)) = 1;
        Offspring = (Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2;
        
        % Polynomial Mutation (applied based on grouping)
        Lower = repmat(Problem.lower,N,1);
        Upper = repmat(Problem.upper,N,1);
        Site = GroupIndex == ChooseGroup;
        mu    = rand(N,1);
        mu    = repmat(mu,1,D);
        temp  = Site & mu<=0.5;
        Offspring       = min(max(Offspring,Lower),Upper);
        Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
            (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
        temp = Site & mu>0.5;
        Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
            (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    end
end