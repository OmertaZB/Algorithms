% OLCA
clc;
clear all;

format long;
format compact;

val_2_reach = 10^(-200);
fhd=@cec17_func; 
runs = 51;
RecordFEsFactor = ...
    [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, ...
    0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
progress = numel(RecordFEsFactor);% 返回14

for problem_size = [10]
    max_nfes = 10000 * problem_size-3000;%最大评价次数
    rand('seed', sum(100 * clock));
    lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];
    fprintf('Running algorithm\n')
    
    %% 随机游走，获取问题特征;预测策略
    H = Rug_inf_entr(problem_size);
    FeatureVector = H;
    predict_label = Ran_For(FeatureVector,problem_size);
    predict_label = str2num(char(predict_label));
    for func = 1 : 30
        %%  使用EDA算法
        if predict_label(func) == 1
            optimum = func * 100.0;
            %% Record the best results
            outcome = [];
            
            fprintf('\n-------------------------------------------------------\n')
            fprintf('Function = %d, Dimension size = %d\n', func, problem_size)
            
            allerrorvals = zeros(progress, runs);%记录14个点
            parfor run_id = 1 : runs
                a = [];
                b = [];
                mu_tm1 = [];
                run_funcvals = [];
                nfes = 0;
                
                %% 参数设置
                PS_max = 100 * problem_size;
                PS_min = 0.5 * (problem_size^2 + problem_size);
                min_pop_size = PS_min;
                max_pop_size = PS_max;
                pop_size = PS_max;
                tao = 0.35;
                nta_max = 5;
                
                popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
                pop = popold;
                fitness = feval(fhd,pop',func);
                fitness = fitness';
                bsf_fit_var = 1e+30;
                bsf_index = 0;
                bsf_solution = zeros(1, problem_size);
                %%%%%%%%%%%%%%%%%%%%%%%% for out
                for i = 1 : pop_size
                    nfes = nfes + 1;
                    if (fitness(i) < bsf_fit_var && isreal(pop(i, :)) && sum(isnan(pop(i, :)))==0 && min(pop(i, :))>=-100 && max(pop(i, :))<=100)
                        bsf_fit_var = fitness(i);
                        bsf_solution = pop(i, :);
                        bsf_index = i;
                    end
                    
                    if nfes > max_nfes;
                        break;
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%% for out
                
                run_funcvals = [run_funcvals;ones(pop_size,1)*bsf_fit_var];
                g = 1;
                while nfes < max_nfes
                    %                     C = {};
                    ct = 0;
                    nta = 0;
                    w = 0;
                    Sup_Pop = [];
                    NN = round(pop_size * tao);
                    %% 选择优势个体（截断选择）
                    [FitnessValue,index]=sort(fitness);
                    i=1:NN;
                    Sup_Pop(i,:) = pop(index(i),:);
                    %% 更新模型 && Algorithm 2                   
                    for i = 1:size( Sup_Pop, 1)
                        a(i,:) = (sum(FitnessValue(1 : NN,1)) - FitnessValue(i)).*  Sup_Pop(i,:);
                        b(i,:) = sum(FitnessValue(1 : NN,1)) - FitnessValue(i);
                    end
                    a = sum(a,1);
                    b = sum(b,1);
                    mut1 = a/b;
                    mut2 = mut1;
                    
                    if g >=2
                        deltat = mut1 - mu_tm1;
                        mut2 = mut1;
                        while (feval(fhd,(mut2 + deltat)',func) < feval(fhd,mut2',func)) && nta < nta_max
                            mut2 = mut2 + nta * deltat;
                            nta = nta + 1;
                        end
                    end
                    mu_tm1 = mut2;
                    
                    ct = (Sup_Pop - repmat(mut2,size(Sup_Pop,1),1))' * (Sup_Pop - repmat(mut2,size(Sup_Pop,1),1));                   
                    ct = ct/size( Sup_Pop, 1);
                    ct = (ct + ct.') / 2;
                    % 建立新的模型（根据以上）
                    pop_new = mvnrnd(mut2,ct,pop_size - 1);
                    pop_new = [pop_new;bsf_solution];
                    pop_new = boundConstraint2(pop_new,popold,lu);
                    %% 计算适应度值
                    fitness3 = feval(fhd,pop_new',func);
                    fitness3 = fitness3';
                    %%%%%%%%%%%%%%%%%%%%%%%% for out
                    for i = 1 : pop_size
                        nfes = nfes + 1;
                        if fitness3(i) < bsf_fit_var 
                            bsf_fit_var = fitness3(i);
                            bsf_solution = pop_new(i, :);
                            bsf_index = i;
                        end
                        
                        if nfes > max_nfes;
                            break;
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%% for out
                    
                    run_funcvals = [run_funcvals;ones(pop_size,1)*bsf_fit_var];
                    
                    fitness = fitness3;
                    pop = pop_new;
                    %% for resizing the population size 调整种群规模
                    
                    %                     plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);
                    %
                    %                     if pop_size > plan_pop_size
                    %                         reduction_ind_num = pop_size - plan_pop_size;
                    %                         if pop_size - reduction_ind_num <  min_pop_size;
                    %                             reduction_ind_num = pop_size - min_pop_size;
                    %                         end
                    %
                    %                         pop_size = pop_size - reduction_ind_num;
                    %                         %SEL = round(ps*pop_size);
                    %                         %fitness_old=fitness;
                    %                         for r = 1 : reduction_ind_num
                    %                             [valBest indBest] = sort(fitness, 'ascend');
                    %                             worst_ind = indBest(end);
                    %                             popold(worst_ind,:) = [];
                    %                             pop(worst_ind,:) = [];
                    %                             fitness(worst_ind,:) = [];
                    %                         end
                    %                     end
                    
                    g = g + 1;
                end% end nfes
                %% Violation Checking
                %             if(max(bsf_solution)>100)
                %                 fprintf('%d th run, Above Max\n', run_id)
                %             end
                %
                %             if(min(bsf_solution)<-100)
                %                 fprintf('%d th run, Below Min\n', run_id)
                %             end
                %
                %             if(~isreal(bsf_solution))
                %                 fprintf('%d th run, Complix\n', run_id)
                %             end
                
                bsf_error_val = abs(bsf_fit_var - optimum);
                %             if bsf_error_val < val_2_reach
                %                 bsf_error_val = 0;
                %             end
                
                if(sum(isnan(bsf_solution))>0)
                    fprintf('%d th run, NaN\n', run_id)
                end
                %% 输出，显示
                fprintf('%d th run, best-so-far error value = %1.8e\n', run_id , bsf_error_val)
                outcome = [outcome bsf_error_val];
                
                %% From Noor Code ( print files )
                errorVals= [];
                for w = 1 : progress
                    bestold = run_funcvals(RecordFEsFactor(w) * max_nfes) - optimum;
                    if abs(bestold)>1e-200
                        errorVals(w)= abs(bestold);
                    else
                        bestold=0;
                        errorVals(w)= bestold;
                    end
                end
                allerrorvals(:, run_id) = errorVals;
            end % end 1 run
            
            %% 使用LSHADE算法
        else
            optimum = func * 100.0;
            %% Record the best results
            outcome = [];
            
            fprintf('\n-------------------------------------------------------\n')
            fprintf('Function = %d, Dimension size = %d\n', func, problem_size)
            
            allerrorvals = zeros(progress, runs);%记录14个点
            parfor run_id = 1 : runs
                %tic
                run_funcvals = [];
                %%  parameter settings
                p_best_rate = 0.11;  
                arc_rate = 1.4;
                memory_size = 5;            
                pop_size = 18 * problem_size;   %18*D   
                goodCR = [];
                goodF = [];
                nfes = 0;
                max_pop_size = pop_size;
                min_pop_size = 4.0;
                ps = 0.5;
                record_CR = [];
                %% 产生初始种群
                popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
                pop = popold;
                fitness = feval(fhd,pop',func);
                fitness = fitness';
                
                bsf_fit_var = 1e+30;
                bsf_index = 0;
                bsf_solution = zeros(1, problem_size);
                
                %%%%%%%%%%%%%%%%%%%%%%%% for out
                for i = 1 : pop_size
                    nfes = nfes + 1;
                    
                    if fitness(i) < bsf_fit_var
                        bsf_fit_var = fitness(i);
                        bsf_solution = pop(i, :);
                        bsf_index = i;
                    end
                    
                    if nfes > max_nfes
                        break;
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%% for out
                run_funcvals = [run_funcvals;ones(pop_size,1)*bsf_fit_var];
                
                %                  memory_sf = 0.3 .* ones(memory_size, 1); %jSO
                %                 memory_cr = 0.8 .* ones(memory_size, 1);  %jSO
                memory_sf = 0.5 .* ones(memory_size, 1);
                memory_cr = 0.5 .* ones(memory_size, 1);
                memory_pos = 1;
                
                archive.NP = arc_rate * pop_size; % the maximum size of the archive
                archive.pop = zeros(0, problem_size); % the solutions stored in te archive
                archive.funvalues = zeros(0, 1); % the function value of the archived solutions
                
                %% main loop  主循环
                
                while nfes < max_nfes
                    
                    pop = popold;
                    [temp_fit, sorted_index] = sort(fitness, 'ascend');
                    
                    mem_rand_index = ceil(memory_size * rand(pop_size, 1));
                    mu_sf = memory_sf(mem_rand_index);
                    mu_cr = memory_cr(mem_rand_index);
                    
                    %% for generating crossover rate   产生交叉率
                    cr = normrnd(mu_cr, 0.1);%正态分布
                    term_pos = find(mu_cr == -1);
                    cr(term_pos) = 0;
                    cr = min(cr, 1);
                    cr = max(cr, 0);
                    
                    %% for generating scaling factor  用于生成缩放因子F
                    sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
                    pos = find(sf <= 0);
                    
                    while ~ isempty(pos)
                        sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
                        pos = find(sf <= 0);
                    end
                    sf = min(sf, 1);
                    r0 = [1 : pop_size];
                    popAll = [pop; archive.pop];
                    [r1, r2] = gnR1R22(pop_size, size(popAll, 1), r0);%随机选两个个体
                    pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions 选择至少两个最好的解决方案
                    randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]    从[1，2，3，...，pNP]中选择
                    randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0 避免rand = 0，ceil（rand）= 0的问题
                    pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions 随机选择前100p％解决方案之一
                    
                    %% 变异
                    
                    %% 变异(不同)
                    vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));   %公式3
                    vi = boundConstraint2(vi, pop, lu);%检查边界
                    
                    %% 交叉
                    mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent 用于指示哪些元素来自父
                    rows = (1 : pop_size)';
                    cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent 选择一个ui的元素不是来自父母的位置
                    jrand = sub2ind([pop_size problem_size], rows, cols);
                    mask(jrand) = false;
                    ui = vi;
                    ui(mask) = pop(mask);
                    
                    children_fitness = feval(fhd, ui', func);
                    children_fitness = children_fitness';
                    
                    bsf_fit_var_old = bsf_fit_var;
                    %%%%%%%%%%%%%%%%%%%%%%%% for out
                    for i = 1 : pop_size
                        nfes = nfes + 1;
                        
                        if children_fitness(i) < bsf_fit_var
                            bsf_fit_var = children_fitness(i);
                            bsf_solution = ui(i, :);
                            bsf_index = i;
                        end
                        
                        if nfes > max_nfes; break; end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%% for out
                    
                    run_funcvals = [run_funcvals;ones(pop_size,1)*bsf_fit_var];
                    dif = abs(fitness - children_fitness);
                    
                    
                    I = (fitness > children_fitness);% 值为1代表父代的解较差，为0代表父代的解较好
                    goodCR = cr(I == 1);
                    goodF = sf(I == 1);
                    dif_val = dif(I == 1);
                    
                    %      isempty(popold(I == 1, :))
                    archive = updateArchive2(archive, popold(I == 1, :), fitness(I == 1));
                    
                    %% 选择 （I == 1: the parent is better; I == 2: the offspring is better）
                    
                    [fitness, I] = min([fitness, children_fitness], [], 2);%选择
                    
                    %run_funcvals = [run_funcvals; fitness];
                    popold(I == 2, :) = ui(I == 2, :);
                    
                    num_success_params = numel(goodCR);
                    
                    if num_success_params > 0
                        sum_dif = sum(dif_val);
                        dif_val = dif_val / sum_dif;
                        
                        %% for updating the memory of scaling factor  用于更新比例因子的存储器
                        memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
                        
                        %% for updating the memory of crossover rate  用于更新交叉速率的内存
                        if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
                            memory_cr(memory_pos)  = -1;
                        else
                            memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
                        end
                        
                        memory_pos = memory_pos + 1;
                        if memory_pos > memory_size;
                            memory_pos = 1;
                        end
                    end
                    %% for resizing the population size 调整种群规模
                    plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);
                    
                    if pop_size > plan_pop_size
                        reduction_ind_num = pop_size - plan_pop_size;
                        if pop_size - reduction_ind_num <  min_pop_size;
                            reduction_ind_num = pop_size - min_pop_size;
                        end
                        
                        pop_size = pop_size - reduction_ind_num;
                        SEL = round(ps*pop_size);
                        fitness_old=fitness;
                        for r = 1 : reduction_ind_num
                            [valBest indBest] = sort(fitness, 'ascend');
                            worst_ind = indBest(end);
                            popold(worst_ind,:) = [];
                            pop(worst_ind,:) = [];
                            fitness(worst_ind,:) = [];
                        end
                        
                        archive.NP = round(arc_rate * pop_size);
                        
                        if size(archive.pop, 1) > archive.NP
                            rndpos = randperm(size(archive.pop, 1));
                            rndpos = rndpos(1 : archive.NP);
                            archive.pop = archive.pop(rndpos, :);
                        end
                    end
                    %                     p_best_rate = ((p_best_rate - pmin)/max_nfes) * nfes + pmin;% jSO
                    %                  g = g + 1;%jSO
                end % nfes
                %% Violation Checking
                if(max(bsf_solution)>100)
                    fprintf('%d th run, Above Max\n', run_id)
                end
                
                if(min(bsf_solution)<-100)
                    fprintf('%d th run, Below Min\n', run_id)
                end
                
                if(~isreal(bsf_solution))
                    fprintf('%d th run, Complix\n', run_id)
                end
                bsf_error_val = abs(bsf_fit_var - optimum);
                %             if bsf_error_val < val_2_reach
                %                 bsf_error_val = 0;
                %             end
                
                if(sum(isnan(bsf_solution))>0)
                    fprintf('%d th run, NaN\n', run_id)
                end
                %% 输出
                fprintf('%d th run, best-so-far error value = %1.8e\n', run_id , bsf_error_val)
                outcome = [outcome bsf_error_val];
                
                %% From Noor Code ( print files )
                errorVals= [];
                for w = 1 : progress
                    bestold = run_funcvals(RecordFEsFactor(w) * max_nfes) - optimum;
                    %                 if abs(bestold)>1e-30
                    errorVals(w)= abs(bestold);
                    %                 else
                    %                     bestold=0;
                    %                     errorVals(w)= bestold;
                    %                 end
                end
                allerrorvals(:, run_id) = errorVals;
                %toc
            end %% end 1 run
        end
        fprintf('\n')
        fprintf('min error value = %1.8e, max = %1.8e, median = %1.8e, mean = %1.2e, std = %1.2e\n', min(outcome), max(outcome), median(outcome), mean(outcome), std(outcome))
    
    end % end 1 func
end
