function [move_points,T] = SparsePointToPlane(source_points,target_points,target_normals,max_icp,max_outer,max_inner,p)
%SPARSE_POINT_TO_PLANE
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
fprintf('sparse point-to-plane:\n');
if length(source_points(:,1))==4
    source_points=source_points(1:3,:);
    target_points=target_points(1:3,:);
end
move_points = source_points;
old_points=source_points;
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
for icp=1:max_icp
    [idx, ~] = knnsearch(NS,move_points','k',1);
    match_normals=target_normals(:,idx);
    match_points= target_points(:,idx);
    fprintf('iteration at %d-%d\n', icp,max_icp);
    for i=1:max_outer
        for j=1:max_inner
            for k=1:length(idx)
                H(k)=match_normals(:,k)'* (move_points(:,k)-match_points(:,k));
            end
            H=H+Lambda/mu;
            % shrink step
            for k=1:length(H(:,1))
                Z(k)=shrink(H(k),p,mu);
            end
            d=Z-Lambda/mu;
            T1=point_to_plane(move_points,match_points,match_normals,d);
            R=SO3(T1);
            t=transl(T1);
            T=T1*T;
            move_points=R*move_points+t';
            dual_max=0;
            for k=1:length(H(:,1))
                dual=norm((old_points(:,k)-move_points(:,k))'*match_normals(:,k));
                if (dual>dual_max)
                    dual_max=dual;
                end
            end
            old_points=move_points;
            if dual_max <1e-5
                break;
            end
        end
        prime_max=0;
        for k=1:length(H(:,1))
            delta(k)=match_normals(:,k)'*(move_points(:,k)-match_points(:,k))-Z(k);
            if (abs(delta(k))>prime_max)
                prime_max=delta(k);
            end
        end
        Lambda=Lambda+mu*delta;
        if dual_max <1e-5 && prime_max<1e-5
            break;
        end
    end
    
end
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
