% generate a 3d coordinate file for a cylindrical calibration object

% prompt = {'Cylinder diameter [mm]','Number of rows:','Number of clumns:','Distance between rows [mm]:','Distance between columns [mm]:'};
% title = 'Input calibration object parameters';
% dims = [1 35];
% % definput = {'20','hsv'};
% answer = inputdlg(prompt,title);
% 
% D=str2num(answer{1});
% Nr=str2num(answer{2});
% Nc=str2num(answer{3});
% dr=str2num(answer{4});
% dc=str2num(answer{5});

D  = 15.3; 
Nr  = 10; 
Nc  = 45; 
dr  = 1; 
dc  = 1; 

% total number of point
Np=Nr*Nc;
% angle between columns
angInc=2*dc/D;

%plot
cFigure; hold all; axisGeom;
P3d=zeros(Nr,Nc,3);
for ir=1:Nr
    z=dr*(ir-1);
    idx_c = 0; 
    for ic=1:Nc
        idx_c = idx_c + 1; 
        ang=angInc*(ic-1)-angInc*(Nc-1)/2-pi/2;
        x=.5*D*cos(ang);
        y=.5*D*sin(ang);
        text(x*(0.95),y*(0.95),z,num2str(ic),'HorizontalAlignment','center',...
            'VerticalAlignment','middle'); 
        P3d(ir,idx_c,:)=[x y z];
    end
end

P3dVec=reshape(P3d,Np,3);
plotV(P3dVec,'sb','MarkerFaceColor','b');

% save file
uisave('P3d','myCylinderCoordinates');


