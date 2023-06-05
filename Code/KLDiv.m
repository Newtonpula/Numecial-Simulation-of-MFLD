function dist=KLDiv(P,Q) % P real, Q empirical

if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end

if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end

%# create an index of the "good" data points
goodIdx = P>0 & Q>0; %# bin counts <0 are not good, either

dist = sum(P(goodIdx) .* log(P(goodIdx) ./Q(goodIdx)));

end
