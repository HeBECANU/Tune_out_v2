function out=std_c4(x)
%https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
%https://www.jstor.org/stable/2682923?origin=crossref&seq=1#metadata_info_tab_contents



x=x(:);
n=numel(x);
if n>1e4
    c4n=1;
else
    c4n=sqrt(2/(n-1))*(gamma(n/2)/gamma((n-1)/2));
    if isnan(c4n)%handle overflow
        %warning('c4 overflow')
        c4n=1;
    end
end

out=std(x)/c4n;
%out=c4n;
end
