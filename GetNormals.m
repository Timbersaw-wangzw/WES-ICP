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