function x = olfx(A,t)
x = max((((A)^t)*[1; zeros(size(A,1)-1,1)]),1);
x = x(1);
end