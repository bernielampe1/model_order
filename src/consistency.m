function [mmdls, maics] = consistency(d, rho, n, beta)
   % generate data with perfect model match
   %X = MGGD_generation(n, d, rho, beta);
   %pd = makedist('beta', 0.5, 0.5);
   %X = random(pd, d, n)-0.5;

   % compute min aic and mdl for each num samples
   mmdls = zeros(n-d, 1);
   maics = zeros(n-d, 1);
   for i = d:n
       X = 100 * rand(d, n)-50;
       [aic, mdl] = model_order(X(:,1:i));
       [~, mmdls(i)] = min(mdl);
       [~, maics(i)] = min(aic);
   end
end
