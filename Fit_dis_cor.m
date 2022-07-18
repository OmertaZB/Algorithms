% clc
% clear all
function B = Fit_dis_cor(problem_size)
fhd = @cec17_func;
% problem_size = 30;
lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];
steps = 3000;
step_size = 10;

% FileName = ['D:\MATLAB2016\bin\bin1\FLDE\CEC2014_BEST\1_10_Best.mat'];
% load(FileName)
% load('2014_1_30_Best.mat');
% load('2015_1_15_Best.mat');
load('2017_1_30_Best.mat');

% b = load('2015_1_15_Best.mat');
% c = load('2017_1_30_Best.mat');
% FF = [a.Optima;b.Optima];
% h = figure(problem_size);
B = [];
% for func = 1 : 15
for func = 1 : 30  
    for T = 1 : 30
        A = 0;
        pop = Random_increasing_walk(problem_size,lu,steps,step_size);
        F = feval(fhd,pop',func);
        F = F';
        for i = 1 : steps
            D(i) = sqrt(sum((pop(i,:) - Optima(1,1:problem_size)).^2));
        end
        
        for j = 1 : steps
            A = A + ((F(j) - mean(F)) * (D(j) - mean(D)));
        end
        Cfd = A/steps;
        FDC(T) = Cfd/(std(F) * std(D));
        FDC = FDC';
    end
    C = mean(FDC);
    B = [B;C]; % B为最终的结果
    
end
%%
% titlename=[num2str(problem_size) 'D'];
% title(titlename,'FontSize',14);
% % figure
% bar(B)
% saveas(h,titlename,'fig');

