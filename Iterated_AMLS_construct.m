function Q = Iterated_AMLS_construct(x,xc,eps,h,mult)
%Function constructs approximation using iterated AMLS approach
Q=zeros(size(x,1),1);
phi = basisFunction_AMLS(x,xc,eps,h);

for i=1:numel(mult)
    Q=Q+phi*mult{i};   
end
end