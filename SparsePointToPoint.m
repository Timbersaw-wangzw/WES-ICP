function [move_points,T,err,rmse] = SparsePointToPoint(source_points,target_points,iter_In,iter_Out,p)
%SPARSE_POINT_TO_POINT 稀疏point to point
%   此处显示详细说明
if length(source_points(:,1))==4
    source_points=source_points(1:3,:);
    target_points=target_points(1:3,:);
end
move_points = source_points;
NS = createns(target_points','NSMethod','kdtree');
lambda_pTpln=zeros(size(move_points));
Z=zeros(size(move_points));
mu=10;
err=[];
rmse=0;
for i=1:iter_Out
    [idx, ~] = knnsearch(NS,move_points','k',1);
    MatchPts= target_points(:,idx);
    for j=1:iter_In
        H=move_points-MatchPts+lambda_pTpln/mu;
        for k=1:length(Z(1,:))
            Z(:,k)=shrink(H(:,k),p,mu);
        end
        C=MatchPts+Z-lambda_pTpln / mu;
        [R,t]=svd_icp(move_points,C);
        move_points=R*move_points+t;
        T=[R,t;[0,0,0,1]];
        T1=eye(4);
        d=norm(T1-T);
        err=[err,d];
        if d<1e-5
            return;
        end
        e=0;
        for k=1:length(move_points(1,:))
            e=e+norm(move_points(:,k)-MatchPts(:,k));
        end
        rmse=sqrt(e/length(move_points(1,:)));
        delta=move_points-MatchPts-Z;
        lambda_pTpln=lambda_pTpln+mu*delta;
    end
end
N=length(move_points(1,:));
RSME=0;
for i=1:N
    RSME=RSME+norm(move_points(:,i)-MatchPts(:,i));
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
