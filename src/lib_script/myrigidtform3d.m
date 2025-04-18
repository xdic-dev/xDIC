function tform = myrigidtform3d(angle,trans)
yaw = angle(1); %z
pitch = angle(2); %y
roll = angle(3); %x
Ryaw = [cosd(yaw),-sind(yaw),0;...
        sind(yaw),cosd(yaw),0;...
        0,       0       ,1]; 
Rpitch = [cosd(pitch), 0,sind(pitch);...
          0,          1,0;...
          -sind(pitch),0,cosd(pitch)]; 
Rroll = [1,0,        0;...
         0,cosd(roll),-sind(roll);...
         0,sind(roll),cosd(roll)];  
R = Ryaw*Rpitch*Rroll; 
tform = rigid3d(R,trans);
end