function [move_points,T,RMSE] = SparseWeightedDistance(source_points,target_points,target_normals,iter_In,iter_Out,p)
%SparseWeightedDistance 稀疏WES-ICP
%   此处显示详细说明
if length(source_points(:,1))==4
    source_points=source_points(1:3,:);
    target_points=target_points(1:3,:);
end
move_points = source_points;
NS = createns(target_points','NSMethod','kdtree');
Lambda=zeros(length(source_points(1,:)),1);
Z=Lambda;
H=Lambda;
delta=Lambda;
mu=10;
T=SE3;
% TargetNormals=GetNormals(target_points,NS);
% load('TargetNormals.mat');
% fileID = fopen('normals.txt','w');
% fprintf(fileID,'%.9f %.9f %.9f\n',TargetNormals);
for i=1:iter_Out
    [idx, ~] = knnsearch(NS,move_points','k',1);
    MatchNormals=target_normals(:,idx);
    MatchPts= target_points(:,idx);
    coeffICP =  1 / (1 + exp(iter_Out - i));
    coeffTDM = 1;
    for j=1:iter_In
        for k=1:length(idx)
            H(k)=MatchNormals(:,k)'* (move_points(:,k)-MatchPts(:,k));
        end
        H=H+Lambda/mu;
        for k=1:length(H(:,1))
            Z(k)=shrink(H(k),p,mu);
        end
        d=Z-Lambda/mu;
        T1=WES_ICP(move_points,MatchPts,MatchNormals,d,coeffICP,coeffTDM);
%         T1=[R,t;[0,0,0,1]];
        R=SO3(T1);
        t=transl(T1);
        T=T1*T;
        move_points=R*move_points+t';
        for k=1:length(H(:,1))
            delta(k)=MatchNormals(:,k)'*(move_points(:,k)-MatchPts(:,k))-Z(k);
        end
        Lambda=Lambda+mu*delta;
    end
end
T=double(T);
N=length(move_points(1,:));
[move_points,T2] = SparsePointToPoint(move_points,target_points,5,20,0.4);
T=T2*T;
RMSE=0;
for i=1:N
    RMSE=RMSE+norm(move_points(:,i)-MatchPts(:,i));
end
RMSE=sqrt(RMSE/N);
end
function T=WES_ICP(SourcePts,TargetPts,TargetNormals,d,w1,w2)
A=zeros(6,6);
b=zeros(6,1);
for i=1:length(SourcePts(1,:))
    hat=dotVec([SourcePts(:,i);1]);
    A1=hat;
    A2=[TargetNormals(:,i);1]'*hat;
    A=A+w1*A1'*A1+w2*(A2'*A2);
    b1=[w1*(SourcePts(:,i)-TargetPts(:,i));0];
    b2=w2*((SourcePts(:,i)-TargetPts(:,i))'*TargetNormals(:,i)-d(i));
    b=b-A1'*b1-A2'*b2;
end
x=A\b;
T=SE3.exp(x);
end
function hat=dotVec(v)
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
function TargetNormals=GetNormals(TargetPts,NS)
TargetNormals=zeros(3,length(TargetPts(1,:)));
for i=1:length(TargetPts(1,:))
    min_neighbors = 10;
    normal = estimateNormal(TargetPts', NS, TargetPts(:,i)', min_neighbors);
    TargetNormals(:,i)=normal;
end
sensorCenter = mean(TargetPts');
for k = 1 : numel(TargetPts)/4
   p1 = sensorCenter(1:3)' - TargetPts(1:3,k);
   p2 = TargetNormals(1:3,k);
   angle = atan2(norm(cross(p1,p2)),p1'*p2);
   if ~(angle > pi/2 || angle < -pi/2)
       TargetNormals(1:3,k)=-1*TargetNormals(1:3,k);
   end
end
end

