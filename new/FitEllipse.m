function [e Cov] = FitEllipse(B, STD);

Mu = mean( B );
X0 = bsxfun(@minus, B, Mu);

%# eigen decomposition [sorted by eigen values]
conf = 2*normcdf(STD)-1;     %# covers around 95% of population
scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions

Cov = cov(X0) * scale;
[V D] = eig(Cov)/length(B);
[D order] = sort(diag(D), 'descend');
D = diag(D);
V = V(:, order);

t = linspace(0,2*pi,100);
e = [cos(t) ; sin(t)];        %# unit circle
VV = V*sqrt(D);               %# scale eigenvectors
e = bsxfun(@plus, VV*e, Mu'); %#' project circle back to orig space