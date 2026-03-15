function value = mat_compare(new,old)

diff = new(:) - old(:);

num = sqrt(sum(abs(diff).^2));
den = sqrt(sum(abs(new(:)).^2));

value = num / den;

end