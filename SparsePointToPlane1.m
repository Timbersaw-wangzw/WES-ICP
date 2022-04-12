 function [move_points,T,RMSE] = SparsePointToPlane(source_points,target_points,target_normals,iter_In,iter_Out,p)
%SPARSE_POINT_TO_POINT 稀疏point to plane
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
old_points=move_points;
% TargetNormals=GetNormals(target_points,NS);
% load('TargetNormals.mat');
% fileID = fopen('normals.txt','w');
% fprintf(fileID,'%.9f %.9f %.9f\n',TargetNormals);
for i=1:iter_Out
    [idx, ~] = knnsearch(NS,move_points','k',1);
    match_bormals=target_normals(:,idx);
    match_points= target_points(:,idx);
    for j=1:iter_In
        for k=1:length(idx)
            H(k)=match_bormals(:,k)'* (move_points(:,k)-match_points(:,k));
        end
        H=H+Lambda/mu;
        for k=1:length(H(:,1))
            Z(k)=shrink(H(k),p,mu);
        end
        d=Z-Lambda/mu;
        for k=1:20
            T1=point_to_plane(move_points,match_points,match_bormals,d);
            R=SO3(T1);
            t=transl(T1);
            T=T1*T;
            move_points=R*move_points+t';
            max=0;
            for m=1:length(H(:,1))
                err=norm(old_points(:,m)-move_points(:,m));
                if (err>max)
                    max=err;
                end
            end
            old_points=move_points;
            if err <1e-5
                break;
            end
        end
    end
    for k=1:length(H(:,1))
        delta(k)=match_bormals(:,k)'*(move_points(:,k)-match_points(:,k))-Z(k);
    end
    Lambda=Lambda+mu*delta;
end
T=double(T);
N=length(move_points(1,:));
RMSE=0;
for i=1:N
    RMSE=RMSE+norm(move_points(:,i)-match_points(:,i));
end
RMSE=sqrt(RMSE/N);
end
function T=point_to_plane(source_points,target_points,target_normals,d)
A=zeros(6,6);
b=zeros(6,1);
for i=1:length(source_points(1,:))
    hat=dotVec([source_points(:,i);1]);
    A1=[target_normals(:,i);1]'*hat;
    A=A+A1'*A1;
    b1=(source_points(:,i)-target_points(:,i))'*target_normals(:,i)-d(i);
    b=b-b1*A1';
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