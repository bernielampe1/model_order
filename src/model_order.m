function [aic, mdl] = model_order(data)

% data: observations vectors assumings each row is a vector
% K: max model order to compute
% return: aic, bic and mdl information criteria vectors for 1:K model order

% compute the correlation matrix
[p, n] = size(data);
R = data * data'./n;

% compute the eigenvectors
e = sort(eig(R), 'descend');

% compute information criteria for all possible k = 1:p-1
aic = zeros(p-1, 1);
mdl = zeros(p-1, 1);
for k = 1:p-1
    aic(k) = gauss_aic(e(k+1:p), k, n, p);
    mdl(k) = gauss_mdl(e(k+1:p), k, n, p);
end

% aic function
function aic = gauss_aic(l, k, n, p)
    aic = -2 * (p-k) * n * log(geomean(l) / mean(l)) + 2 * k * (2 * p - k);
end

% mdl function
function mdl = gauss_mdl(l, k, n, p)
    mdl = -(p-k) * n * log(geomean(l) / mean(l)) + 0.5 * k * (2 * p - k) * log(n);
end

end
