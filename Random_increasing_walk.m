% Random increasing walk algorithm
function walk = Random_increasing_walk(D,lu,steps,step_size)
% problem_size = 10;
% D = problem_size;
% lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];
% steps = 3000;
% step_size = 10;
walk = zeros(steps,D);
count = 1;

start_point = repmat(lu(1, :), 1, 1) + rand(1, D) .* (repmat(lu(2, :) - lu(1, :), 1, 1));
walk(1,:) = start_point;

while count < steps
    for i = 1 : D
        bb = 10 * rand;
        walk(count+1,i) = walk(count,i) + bb;
        if walk(count + 1,i) > 100
            walk(count+1,i) = walk(count+1,i) - 200;
        end
    end
    count = count + 1;
end