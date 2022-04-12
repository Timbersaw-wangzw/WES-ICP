clc
clear
addpath('github_repo');
data01=importdata('coarse_source_points0.6.txt');
data02=importdata('target_points.txt');
source_points=data01';
target_points_normals=data02';
target_points=target_points_normals(1:3,:);
target_normals=target_points_normals(4:6,:);

% load('source_points.mat')
% load('target_points.mat')
% load('target_normals.mat')

hold on
% plot3(SourcePts(1,:),SourcePts(2,:),SourcePts(3,:),'r.');
% plot3(TargetPts(1,:),TargetPts(2,:),TargetPts(3,:),'b.');
axis off
% plot3(source_points(1,:),source_points(2,:),source_points(3,:),'b.','MarkerSize',2.5);
% [move_points,q,~,e2]=SparsePointToPoint(source_points,target_points,5,100,0.4);
[move_points,T,RMSE]=SparsePointToPlane(source_points,target_points,target_normals,50,20,5,0.4);
% [move_points,T,RMSE]=SparseWeightedDistance(source_points,target_points,target_normals,50,100,5,0.4);
% hold on

% plot3(target_points(1,:),target_points(2,:),target_points(3,:),'r.','MarkerSize',2.5);
% plot3(move_points(1,:),move_points(2,:),move_points(3,:),'b.','MarkerSize',2.5);
plot3(target_points(1,:),target_points(2,:),target_points(3,:),'r.');
plot3(move_points(1,:),move_points(2,:),move_points(3,:),'b.');
%% bunny registration
% BunnyICP;
%% highspeed train
% HighSpeedTrainICP;
%% blade surface1
% BladeICP;
% [MovePts,e1]=IRLS_ICP(SourcePts(1:3,:),TargetPts(1:3,:),0.0001,'Cauchy',200);
% [MovePts,T] = icp(SourcePts,TargetPts,100);
% [MovePts,q,~,e2]=SparsePointToPoint(SourcePts,TargetPts,5,100,0.4);
% [MovePts,q,e]=SparsePointToPlane(SourcePts,TargetPts,5,100,0.4);
% [MovePts,T,err,rmse] = sparse_point_to_point(SourcePts(1:3,:),TargetPts(1:3,:),5,100,0.3);
% R=[0.87302 0.0742131 -0.482005;
% -0.261103  0.905894 -0.333438;
%    0.4119  0.416951  0.810241;];
% t=[-0.0392534;-0.0115146;-0.0697322;];
% MovePts=R*SourcePts+t;
% R=[ 0.951495  0.127856 -0.279838;
% -0.168169  0.977794 -0.125055;
%  0.257635  0.166049  0.951868;];
% t=[-0.0660745;-0.0220316;-0.162078;];
% MovePts=R*SourcePts+t;
% hold on
% axis off
% plot3(SourcePts(1,:),SourcePts(2,:),SourcePts(3,:),'b.','MarkerSize',2.5);
% plot3(TargetPts(1,:),TargetPts(2,:),TargetPts(3,:),'r.','MarkerSize',2.5);
% plot3(MovePts(1,:),MovePts(2,:),MovePts(3,:),'b.','MarkerSize',2.5);
