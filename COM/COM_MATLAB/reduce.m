function result = reduce(var1)
% --- Reduce 1x1xn array to 1xn (aka squeeze)
out = zeros(1,length(var1));
out(1,:) = var1(1,1,:);
result=out;