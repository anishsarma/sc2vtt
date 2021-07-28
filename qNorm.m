function y = qNorm(x,p)
y = zeros(1,size(x,2));
for tt = 1:size(x,2)

    y(tt) = norm(x(:,1:tt),p);

end