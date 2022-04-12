function M = HomoOdot(p)
% operate 'odot' homogenous point p into a 4x6 matrix 
M=zeros(4,6);
M(1:3,1:3)=eye(3);
M(1:3,4:6)=-1*skew(p(1:3));
end

