clc
clear
addpath('github_repo');
data01=importdata('coarse_source_points0.5.txt');
data02=importdata('target_points.txt');
source_points=data01'; 
target_points_normals=data02';
target_points=target_points_normals(1:3,:);
target_normals=target_points_normals(4:6,:);

% 
% [move_points,T]=SparsePointToPoint(sour ce_points,target_points,50,20,5,0.4);
% [move_points,T]=SparsePointToPlane(source_points,target_points,target_normals,50,20,5,0.4);
[move_points,T]=SparseWeightedDistance(source_points,target_points,target_normals,50,20,5,0.4);

hold on
axis off
% plot3(source_points(1,:),source_points(2,:),source_points(3,:),'b.','MarkerSize',2.5);
plot3(target_points(1,:),target_points(2,:),target_points(3,:),'r.','MarkerSize',2.5);
plot3(move_points(1,:),move_points(2,:),move_points(3,:),'b.','MarkerSize',2.5);
