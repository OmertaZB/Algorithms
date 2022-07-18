% clc
% clear all
function HHE =  Rug_inf_entr(problem_size)
fhd = @cec17_func;
% fhd = @cec15_func;

% problem_size = 30;
lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];
steps = 3000;
step_size = 10;

% load('2014_1_30_Best.mat');
% load('2015_1_15_Best.mat');
load('2017_1_30_Best.mat');


% for func = 1 : 15
HHE = [];
for func = 1 : 30
    
    G = 1;
    for T = 1 : 30
        %         A = 0;
        pop = Random_increasing_walk(problem_size,lu,steps,step_size);
        F = feval(fhd,pop',func);
        F = F';
        
        for i = 2 : steps
            e(i-1) = (F(i) - F(i-1));
        end
        e = max(e);
        
        g = 1;
        
        for e = [0 e/128 e/64 e/32 e/16 e/8 e/4 e/2 e]
            count_m1_0 = 0;
            count_m1_1 = 0;
            count_0_1 = 0;
            count_0_m1 = 0;
            count_1_m1 = 0;
            count_1_0 = 0;
            for i = 2 : steps
                if (F(i) - F(i-1)) < -e
                    Se(i) = -1;
                elseif abs(F(i) - F(i-1)) <= e
                    Se(i) = 0;
                elseif (F(i) - F(i-1)) > e
                    Se(i) = 1;
                end
            end
            Se(1) = [];
            for i = 1 : steps - 2
                if (Se(i) == -1) && (Se(i + 1) == 0)
                    count_m1_0 = count_m1_0 + 1;
                elseif (Se(i) == -1) && (Se(i + 1) == 1)
                    count_m1_1 = count_m1_1 + 1;
                elseif (Se(i) == 0) && (Se(i + 1) == 1)
                    count_0_1 = count_0_1 + 1;
                elseif (Se(i) == 0) && (Se(i + 1) == -1)
                    count_0_m1 = count_0_m1 + 1;
                elseif (Se(i) == 1) && (Se(i + 1) == -1)
                    count_1_m1 = count_1_m1 + 1;
                elseif (Se(i) == 1) && (Se(i + 1) == 0)
                    count_1_0 = count_1_0 + 1;
                end
            end
            p_count_m1_0 = count_m1_0/steps;
            p_count_m1_1 = count_m1_1/steps;
            p_count_0_1 = count_0_1/steps;
            p_count_0_m1 = count_0_m1/steps;
            p_count_1_m1 = count_1_m1/steps;
            p_count_1_0 = count_1_0/steps;
            
            
            
            if p_count_m1_0 == 0
                He1 = 0;
            else
                He1 = p_count_m1_0 * (log(p_count_m1_0)/log(6));
            end
            
            if p_count_m1_1 == 0
                He2 = 0;
            else
                He2 = p_count_m1_1 * (log(p_count_m1_1)/log(6));
            end
            if p_count_0_1 == 0
                He3 = 0;
            else
                He3 = p_count_0_1 * (log(p_count_0_1)/log(6));
            end
            if p_count_0_m1 == 0
                He4 = 0;
            else
                He4 = p_count_0_m1 * (log(p_count_0_m1)/log(6));
            end
            if p_count_1_m1 == 0
                He5 = 0;
            else
                He5 = p_count_1_m1 * (log(p_count_1_m1)/log(6));
            end
            if p_count_1_0 == 0
                He6 = 0;
            else
                He6 = p_count_1_0 * (log(p_count_1_0)/log(6));
            end
            HE(G,g) = - sum(He1 + He2 + He3 + He4 + He5 + He6);
%             HE(find(HE == inf)) = 0;
            g = g + 1;
        end
        
        G =G + 1;
    end
    HHE = [HHE;mean(HE,1)];
end

%%
% hold on
% plot([0 1 2 3 4 5 6 7 8],HHE)
