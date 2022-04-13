function [move_points,T] = SparsePointToPoint(source_points,target_points,max_icp,max_outer,max_inner,p)
%SPARSE_POINT_TO_POINT sparse point to point
%   INPUT:
%   source_points: source point clouds
%   target_points: target point clouds
%   max_icp: maximum iteration number of searching cloest points and
%   solving R and t
%   max_out: maximum iteration number of solving ALM by ADMM
%   max_in: maximum iteration number of solving R and t
%   p: p-norm
%   OUTPUT:
%   move_points: final results
%   T: optimal homogenous transformation matrix 
fprintf('sparse point-to-point:\n');
if length(source_points(:,1))==4
    source_points=source_points(1:3,:);
    target_points=target_points(1:3,:);
end
move_points = source_points;
NS = createns(target_points','NSMethod','kdtree');
Num=length(source_points(1,:));
old_points=source_points;
lambda=zeros(size(move_points));
Z=zeros(size(move_points));
mu=10;
T=eye(4);
for icp=1:max_icp
    [idx, ~] = knnsearch(NS,move_points','k',1);
    match_points= target_points(:,idx);
    fprintf('iteration at %d-%d\n', icp,max_icp);
    for i=1:max_outer
        for j=1:max_inner
            H=move_points-match_points+lambda/mu;
            for k=1:Num
                Z(:,k)=shrink(H(:,k),p,mu);
            end
            C=match_points+Z-lambda / mu;
            [R,t]=svd_icp(move_points,C);
            T1=eye(4);
            T1(1:3,1:3)=R;
            T1(1:3,4)=t;
            T=T1*T;
            move_points=R*move_points+t;
            dual_max=0;
            for m=1:Num
                dual=norm((old_points(:,m)-move_points(:,m)));
                if (dual>dual_max)
                    dual_max=dual;
                end
            end
            old_points=move_points;
            if dual_max <1e-5
                break;
            end
        end
        delta=move_points-match_points-Z;
        prime_max=0;
        for k=1:Num
            prime=norm(move_points(:,k)-match_points(:,k)-Z(:,k));
            if (prime>prime_max)
                prime_max=prime;
            end
        end
        lambda=lambda+mu*delta;
        if dual_max <1e-5 && prime_max<1e-5
            break;
        end
    end
end
end
function [R,t]=svd_icp(source_points,target_points)
centerA=mean(source_points')';
centerB=mean(target_points')';
tempA=source_points-centerA;
tempB=target_points-centerB;
H=zeros(3,3);
for i=1:length(source_points(1,:))
    H=H+tempA(:,i)*tempB(:,i)';
end
[U,~,V]=svd(H);
R=V*U';
t=centerB-R*centerA;
end
