classdef DSGEA < ALGORITHM
% <multi> <real/binary> <large/none> <constrained/none> <sparse>
% Dynamic Sparse Grouping Evolutionary Algorithm

%------------------------------- Reference --------------------------------
% Z. Yingjie, Y. Liu, J. Zou, S. Yang, and J. Zheng. "An evolutionary
% algorithm based on dynamic sparse grouping for sparse large scale
% multiobjective optimization." Information Sciences 631 (2023): 449-467.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem) 
            %% Parameter setting
            TDec    = [];
            TMask   = [];
            TempPop = [];
            Fitness = zeros(1,Problem.D);
            REAL    = ~strcmp(Problem.encoding,'binary');
            
            %% Initialization for fitness evaluation
            for i = 1 : 1+4*REAL            
                if REAL
                    Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                else                    
                    Dec = ones(Problem.D,Problem.D);
                end
                Mask       = eye(Problem.D);
                Population = Problem.Evaluation(Dec.*Mask);
                TDec       = [TDec;Dec];
                TMask      = [TMask;Mask];
                TempPop    = [TempPop,Population];
                Fitness    = [Fitness;NDSort([Population.objs,Population.cons],inf)];                
            end
            Fitness = median(Fitness);
            
            %% Generate initial population
            if REAL
                Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));                
            else
                Dec = ones(Problem.N,Problem.D);
            end
            Mask = zeros(Problem.N,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;                
            end
            
            Population = Problem.Evaluation(Dec.*Mask);
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,TempPop],[Dec;TDec],[Mask;TMask],Problem.N);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                iter = Problem.FE/Problem.maxFE;
                
                % Algorithm 3: Fitness Update
                Fitness = FitnessUpdate(Fitness,Population,iter);
                
                % Mating pool selection
                MatingPool = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);   
                
                % Algorithm 5: Variation
                [OffDec,OffMask] = Variation(Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness,REAL,iter);
                
                % Evaluation and environmental selection
                Offspring = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
                            
            end
        end
    end
end