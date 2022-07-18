function vi = boundConstraint2 (vi, pop, lu)

% if the boundary constraint is violated, set the value to be the middle
% of the previous value and the bound


[NP, ~] = size(pop);  % the population size and the problem's dimension

%% check the lower bound
xl = repmat(lu(1, :), NP, 1);
pos = vi < xl;
vi(pos) = (pop(pos) + xl(pos)) / 2;% ���в���

%% check the upper bound
xu = repmat(lu(2, :), NP, 1);
pos = vi > xu;
vi(pos) = (pop(pos) + xu(pos)) / 2;