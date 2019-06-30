function out=frac_diff(x,y,meth)
%calculate the fractional difference
% TODO
% document
% test
% more methods eg take min or max as denominator
out=(x-y)./((x+y)/2);
end