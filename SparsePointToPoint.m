function [MovePts,T,err,rmse] = SparsePointToPoint(source_points,target_points,iter_In,iter_Out,p)
%SPARSE_POINT_TO_POINT 稀疏point to point
%   此处显示详细说明
if length(source_points(:,1))==4
    source_points=source_points(1:3,:);
    target_points=target_points(1:3,:);
end
MovePts = source_points;
NS = createns(target_points','NSMethod','kdtree');
Lambda=zeros(size(MovePts));
mu=10;
err=[];
rmse=0;
for i=1:iter_Out
    [idx, ~] = knnsearch(NS,MovePts','k',1);
    MatchPts= target_points(:,idx);
    for j=1:iter_In
        H=MovePts-MatchPts+Lambda/mu;
        Z=zeros(size(MovePts));
        for k=1:length(Z(1,:))
            Z(:,k)=shrink(H(:,k),p,mu);
        end
        C=MatchPts+Z-Lambda / mu;
        [R,t]=svd_icp(MovePts,C);
        MovePts=R*MovePts+t;
        T=[R,t;[0,0,0,1]];
        T1=eye(4);
        d=norm(T1-T);
        err=[err,d];
        if d<1e-5
            return;
        end
        e=0;
        for k=1:length(MovePts(1,:))
            e=e+norm(MovePts(:,k)-MatchPts(:,k));
        end
        rmse=sqrt(e/length(MovePts(1,:)));
        delta=MovePts-MatchPts-Z;
        Lambda=Lambda+mu*delta;
    end
end
N=length(MovePts(1,:));
RSME=0;
for i=1:N
    RSME=RSME+norm(MovePts(:,i)-MatchPts(:,i));
end
RSME=sqrt(RSME/N);
end
function [R,t]=svd_icp(SourcePts,TargetPts)
centerA=mean(SourcePts')';
centerB=mean(TargetPts')';     
tempA=SourcePts-centerA;        
tempB=TargetPts-centerB;
H=zeros(3,3);
for i=1:length(SourcePts(1,:))
    H=H+tempA(:,i)*tempB(:,i)';
end
[U,~,V]=svd(H);
R=V*U';
t=centerB-R*centerA;
end
