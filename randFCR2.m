function [F,CR] = randFCR2(NP, CRm, CRsigma, Fm,  Fsigma)       %NP = 100, CRm = 0.5, Fm = 0.1, Fsigma = 0.1
%% generate CR
CR = CRm + CRsigma*randn(NP, 1);
CR = min(1, max(0, CR));                % truncated to [0 1]
%% generate F
F = randCauchy(NP, 1, Fm, Fsigma);          
F = min(1, F);
% we don't want F = 0. So, if F<=0, we regenerate F (instead of trucating it to 0)
pos = find(F <= 0);
while ~ isempty(pos)
    F(pos) = randCauchy(length(pos), 1, Fm, Fsigma);
    F = min(1, F);                      % truncation
    pos = find(F <= 0);
end
% Cauchy distribution: cauchypdf = @(x, mu, delta) 1/pi*delta./((x-mu).^2+delta^2)
%% randCauchy产生F的函数
function result = randCauchy(m, n, mu, delta)
result = mu + delta * tan(pi * (rand(m, n) - 0.5));





