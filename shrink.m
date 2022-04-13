function out = shrink(h,p,mu)
%SHRINK shrink algorithm
alpha_a=((2. / mu)*(1. - p))^(1. / (2. - p));
hTilde= alpha_a + (p / mu)*alpha_a^(p - 1);
hNorm=norm(h);
if hNorm<=hTilde
    out=0;
    return 
else
    beta=((alpha_a ) / hNorm + 1.) / 2.;
    for i=1:3
        beta = 1-(p / mu)*(hNorm^(p - 2))*(beta^(p - 1));
    end 
    out= beta*h;
    return 
end

