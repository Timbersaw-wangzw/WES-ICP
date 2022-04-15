function [move_points,T] = SparseWeightedDistance(source_points,target_points,target_normals,max_icp,max_outer,max_inner,p)
%SparseWeightedDistance WES-ICP
%   INPUT:
%   source_points: source point clouds
%   target_points: target point clouds
%   target_normals: target point clouds normals
%   max_icp: maximum iteration number of searching cloest points and
%   solving R and t
%   max_out: maximum iteration number of solving ALM by ADMM
%   max_in: maximum iteration number of solving R and t
%   p: p-norm
%   OUTPUT:
%   move_points: final results
%   T: optimal homogenous transformation matrix 
fprintf('WES-ICP:\n');
if length(source_points(:,1))==4
    source_points=source_points(1:3,:);
    target_points=target_points(1:3,:);
end
move_points = source_points;
NS = createns(target_points','NSMethod','kdtree');
Num=length(source_points(1,:));
old_points=source_points;
%_pTpln denotes the point to plane
%_pTp denotes the point to point
lambda_pTpln=zeros(Num,1);
Z_pTpln=lambda_pTpln;
H_pTpln=lambda_pTpln;
delta_pTpln=lambda_pTpln;
lambda_pTp=zeros(size(move_points));
Z_pTp=zeros(size(move_points));
mu=10;
T=SE3;
for icp=1:max_icp
    [idx, ~] = knnsearch(NS,move_points','k',1);
    match_normals=target_normals(:,idx);
    match_points= target_points(:,idx);
    % weight the point-to-point and point-to-plane distance
    w_pTp =  2 / (1+exp(max_icp - icp));
    w_pTpln = 1-w_pTp;
%     if icp==max_icp-10
%         w_pTp =  1;
%         w_pTpln = 0;
%     end
    fprintf('iteration at %d-%d\n', icp,max_icp);
    for i=1:max_outer
        for j=1:max_inner
            % sparse point-to-point shrink
            H_point_to_point=move_points-match_points+lambda_pTp/mu;
            for  k=1:Num
                Z_pTp(:,k)=shrink(H_point_to_point(:,k),p,mu);
            end
            inter_points=match_points+Z_pTp-lambda_pTp / mu;
            % sparse point-to-plane shrink
            for  k=1:Num
                H_pTpln(k)=match_normals(:,k)'* (move_points(:,k)-match_points(:,k));
            end
            H_pTpln=H_pTpln+lambda_pTpln/mu;
            for  k=1:Num
                Z_pTpln(k)=shrink(H_pTpln(k),p,mu);
            end
            inter_d=Z_pTpln-lambda_pTpln/mu;
            % solve R and t
            T1=WES_ICP(move_points,inter_points,match_points,match_normals,inter_d,w_pTp,w_pTpln);
            R=SO3(T1);
            t=transl(T1);
            T=T1*T;
            move_points=R*move_points+t';
            dual_max=0;
            for  k=1:Num
                dual_pTpln=norm((old_points(:,k)-move_points(:,k))'*match_normals(:,k));
                dual_pTp=norm((old_points(:,k)-move_points(:,k)));
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
        delta_pTp=move_points-match_points-Z_pTp;
        prime_max=0;
        for  k=1:Num
            delta_pTpln(k)=match_normals(:,k)'*(move_points(:,k)-match_points(:,k))-Z_pTpln(k);
            prime_pTpln=abs(delta_pTpln(k));
            prime_pTp=norm(move_points(:,k)-match_points(:,k)-Z_pTp(:,k));
            prime=max(prime_pTpln,prime_pTp);
            if (prime>prime_max)
                prime_max=prime;
            end
        end
        lambda_pTpln=lambda_pTpln+mu*delta_pTpln;
        lambda_pTp=lambda_pTp+mu*delta_pTp;
        if dual_max <1e-5 && prime_max<1e-5
            break;
        end
    end
end
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


