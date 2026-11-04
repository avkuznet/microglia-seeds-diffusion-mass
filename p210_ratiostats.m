function R = p210_ratiostats(num, den, modelRatio, nboot, nMice, rho)
%P210_RATIOSTATS  Compare a model ratio with a measured ratio of two samples.
%
%   R = P210_RATIOSTATS(num, den, modelRatio)
%   R = P210_RATIOSTATS(num, den, modelRatio, nboot)
%
%   num and den are the individual-animal values at the two time points or
%   conditions; they need not be paired or of equal length.
%
%   WHY MORE THAN ONE STATISTIC.  Graft amyloid loads are strongly
%   right-skewed - in the Source Data the individual grafts span a factor of
%   fifty - so propagating standard errors through a ratio of means, which
%   assumes normality, understates the lower tail badly and can turn a
%   marginal discrepancy into an apparently decisive one.  The logarithmic
%   scale, on which these data are approximately normal, is the appropriate
%   one; the bootstrap on the ratio of means is reported alongside it as a
%   check that does not assume a distribution at all.
%
%   nMice = [n1 n2] gives the number of ANIMALS behind the two samples, when
%   the measurements are grafts and a mouse can contribute more than one.
%   d'Errico et al. report every graft-level panel as "one graft from n
%   mice", so the values are clustered and the naive standard error is too
%   small.  With mean cluster size mbar the design effect is 1+(mbar-1)rho;
%   rho defaults to 1, the worst case, which bounds the inflation rather
%   than estimating it.  Passing nMice leaves the point estimates alone and
%   widens the log-scale interval only.
%
%   Returns R with fields: rMean, seMeanGauss, sigGauss, rMedian,
%   logRatio, logSE, logT, logCI, sigLog, bootCI, bootP, nboot, deff.
%   The log-scale fields are NaN if either sample contains a zero.

if nargin < 4 || isempty(nboot), nboot = 40000; end
if nargin < 5, nMice = []; end
if nargin < 6 || isempty(rho), rho = 1; end
num = num(:); den = den(:);
n1 = numel(num); n2 = numel(den);

% ---- ratio of means, with Gaussian propagation (reported, not trusted) --
m1 = mean(num);  m2 = mean(den);
e1 = std(num)/sqrt(n1);  e2 = std(den)/sqrt(n2);
R.rMean       = m1/m2;
R.seMeanGauss = R.rMean * hypot(e1/m1, e2/m2);
R.sigGauss    = (modelRatio - R.rMean) / R.seMeanGauss;

% ---- ratio of medians ---------------------------------------------------
R.rMedian = median(num) / median(den);

% ---- logarithmic scale --------------------------------------------------
if all(num > 0) && all(den > 0)
    l1 = log(num);  l2 = log(den);
    R.logRatio = mean(l1) - mean(l2);
    if isempty(nMice)
        R.deff = [1 1];
    else
        R.deff = 1 + ([n1 n2]./nMice(:)' - 1) * rho;   % Kish design effect
    end
    R.logSE    = sqrt(var(l1)/n1*R.deff(1) + var(l2)/n2*R.deff(2));
    R.logT     = R.logRatio / R.logSE;                 % is the change itself real?
    R.logCI    = exp(R.logRatio + [-1.96 1.96]*R.logSE);
    R.sigLog   = (log(modelRatio) - R.logRatio) / R.logSE;
else
    [R.logRatio, R.logSE, R.logT, R.sigLog] = deal(NaN);
    R.logCI = [NaN NaN];  R.deff = [NaN NaN];
end

% ---- bootstrap on the ratio of means ------------------------------------
s = rng; rng(7);                      % fixed seed: the printed numbers repeat
b = zeros(nboot,1);
for k = 1:nboot
    b(k) = mean(num(randi(n1,n1,1))) / mean(den(randi(n2,n2,1)));
end
rng(s);
b = sort(b);
R.bootCI = [b(max(1,floor(0.025*nboot))) b(min(nboot,ceil(0.975*nboot)))];
R.bootP  = mean(b >= modelRatio);     % one-sided, against the model value
R.nboot  = nboot;
end
