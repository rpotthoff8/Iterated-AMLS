%For multi Dimensions
function [phi,h] = basisFunction_AMLS(x,xc,eps,h)
%From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer
%and Zhang, Equation 22

s=size(x,2); %Dimensions
dist=zeros(size(xc,1),size(x,1),s);
for i=1:s
    dist(:,:,i)=x(:,i)'-xc(:,i); %solve distance in each dimension, Eqn 16 from Javier's notes
end
R_square=sum(dist.^2,3); %Eqn 17 from Javier's notes

if h==0
    R_square_find=R_square;
    R_square_find(R_square_find<=0)=nan;
    h=sqrt(max(min(R_square_find)));
end

phi=((eps^s)/(sqrt(pi^s)))*exp(-((eps^2)*(R_square'))/h^2); %Eqn 10 from 'Development of RBF-DQ method... Y.L Wu

end
