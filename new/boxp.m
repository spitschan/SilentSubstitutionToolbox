function boxp(x);
n = length(x);

% First get the quartiles.
q = quantile(x, [0.25 0.5 0.75]);
% Find the interquartile range.
iq = q(3) - q(1);
% Find the outer limits.
UL = q(3) + 1.5*iq;
LL = q(1) - 1.5*iq;

% Find any outliers.
ind = [find(x > UL) find(x < LL)];
outs = x(ind);
% Get the adjacent values. Find the
% points that are NOT outliers.
inds = setdiff(1:n,ind);
% Get their min and max.
adv = [x(inds(1)) x(inds(end))];
    
% Now draw the necessary pieces.
% Draw the quartiles.
plot([1 3],[q(1),q(1)], '-k')
hold on
plot([1 3],[q(2),q(2)], '-k')
plot([1 3],[q(3),q(3)], '-k')
% Draw the sides of the box
plot([1 1],[q(1),q(3)], '-k')
plot([3 3],[q(1),q(3)], '-k')
% Draw the whiskers.
plot([2 2],[q(1),adv(1)], '-k',[1.75 2.25],[adv(1) adv(1)], '-k')
plot([2 2],[q(3),adv(2)], '-k', [1.75 2.25],[adv(2) adv(2)], '-k')
% Now draw the outliers with symbols.
plot(2*ones(size(outs)), outs,'.k');

g(1,1)=gramm('x',cars.Origin_Region,'y',cars.Horsepower,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5);
g(1,2)=copy(g(1));
g(2,1)=copy(g(1));
g(2,2)=copy(g(1));