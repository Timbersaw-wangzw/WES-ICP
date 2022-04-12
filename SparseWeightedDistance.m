function [move_points,T,RMSE] = SparseWeightedDistance(source_points,target_points,target_normals,iter_In,iter_Out,p)
%SparseWeightedDistance 稀疏WES-ICP
%   此处显示详细说明
if length(source_points(:,1))==4
    source_points=source_points(1:3,:);
    target_points=target_points(1:3,:);
end
move_points = source_points;
NS = createns(target_points','NSMethod','kdtree');

lambda_pTpln=zeros(length(source_points(1,:)),1);
Z_pTpln=lambda_pTpln;
H_pTpln=lambda_pTpln;
delta_pTpln=lambda_pTpln;

lambda_pTp=zeros(size(move_points));
Z_pTp=zeros(size(move_points));
mu=10;
T=SE3;
for i=1:iter_Out
    [idx, ~] = knnsearch(NS,move_points','k',1);
    match_normals=target_normals(:,idx);
    match_points= target_points(:,idx);
    w_ICP =  1 / (1 + exp(iter_Out - i));
    w_TDM = 1;
    for j=1:iter_In
       % sparse point-to-point shrink
        H_point_to_point=move_points-match_points+lambda_pTp/mu;
        for k=1:length(Z_pTp(1,:))
            Z_pTp(:,k)=shrink(H_point_to_point(:,k),p,mu);
        end
        inter_points=match_points+Z_pTp-lambda_pTp / mu;
       % sparse point-to-plane shrink
        for k=1:length(idx)
            H_pTpln(k)=match_normals(:,k)'* (move_points(:,k)-match_points(:,k));
        end
        H_pTpln=H_pTpln+lambda_pTpln/mu;
        for k=1:length(H_pTpln(:,1))
            Z_pTpln(k)=shrink(H_pTpln(k),p,mu);
        end
        
        inter_d=Z_pTpln-lambda_pTpln/mu;
        T1=WES_ICP(move_points,inter_points,match_points,match_normals,inter_d,w_ICP,w_TDM);
        R=SO3(T1);
        t=transl(T1);
        T=T1*T;
        move_points=R*move_points+t';
        dual_max=0;
        for m=1:length(H_pTpln(:,1))
            dual_pTpln=norm((old_points(:,m)-move_points(:,m))'*match_normals(:,m));
            dual_pTp=norm((old_points(:,m)-move_points(:,m)));
            dual=max(dual_pTpln,dual_pTp);
            if (dual>dual_max)
                dual_max=dual;
            end
        end
        old_points=move_points;
        if dual_max <1e-5
            break;
        end
    end
    for k=1:length(H_pTpln(:,1))
            delta_pTpln(k)=match_normals(:,k)'*(move_points(:,k)-match_points(:,k))-Z_pTpln(k);
    end
    lambda_pTpln=lambda_pTpln+mu*delta_pTpln;
end
T=double(T);
N=length(move_points(1,:));
[move_points,T2] = SparsePointToPoint(move_points,target_points,5,20,0.4);
T=T2*T;
RMSE=0;
for i=1:N
    RMSE=RMSE+norm(move_points(:,i)-match_points(:,i));
end
RMSE=sqrt(RMSE/N);
end
function T=WES_ICP(source_points,inter_points,target_points,target_normals,d,w1,w2)
A=zeros(6,6);
b=zeros(6,1);
for i=1:length(source_points(1,:))
    hat=dotVec([source_points(:,i);1]);
    A1=hat;
    A2=[target_normals(:,i);0]'*hat;
    A=A+w1*(A1'*A1)+w2*(A2'*A2);
    b1=[w1*(source_points(:,i)-inter_points(:,i));0];
    b2=w2*((source_points(:,i)-target_points(:,i))'*target_normals(:,i)-d(i));
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

