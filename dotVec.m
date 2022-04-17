function hat=dotVec(v)
% INPUT:
%   v: 4x1 vector
% OUTPUT:
%   hat:4x6 matrix
hat=zeros(4,6);
hat(1:3,1:3)=eye(3);
sysm=zeros(3,3);
sysm(1,2)=-1*v(3);
sysm(1,3)=v(2);
sysm(2,1)=v(3);
sysm(2,3)=-1*v(1);
sysm(3,1)=-1*v(2);
sysm(3,2)=v(1);
hat(1:3,4:6)=-1*sysm;
end