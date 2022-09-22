function [ res ] = Ainverse( par,f,DtD,tau_n)

A = (par.ep*tau_n)*DtD + 1;
g=fftn(f);
res=real(ifftn(g./A));
end

