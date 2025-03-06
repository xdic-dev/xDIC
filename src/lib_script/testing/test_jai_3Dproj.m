frame_target = 150; 
mm=2; 
z_jai = zeros(length(flxyfilt),Nframe); 

xyJAI = flip(((flxyfilt-flip(c/2))/res),2); 
xyJAI(:,1,:) = xyJAI(:,1,:); 
[xq,yq] = meshgrid(-10:0.05:10,-8:0.05:18);

for iframe = 1:150 %:Ntrial 
   PCcurrent = T{mm}.FC.xyz(:,:,iframe); 
   validLogic = ~isnan(PCcurrent(:,1));
   % test interpolation 
   F = scatteredInterpolant(PCcurrent(validLogic,1),PCcurrent(validLogic,2),PCcurrent(validLogic,3));
   F.Method = 'natural';
   vq1 = F(xq,yq);
   
   xyinterpolated = [xq(:),yq(:)]; 
   k = dsearchn(xyinterpolated,xyJAI(:,:,iframe)); 
   
   z_jai(:,iframe) = vq1(k);    
end
%%

size(xyJAI)
size(z_jai)
z_jai_format = zeros(length(xyJAI),1,Nframe);
z_jai_format(:,1,:) = z_jai;
p = cat(2,xyJAI,z_jai_format);
pcell = squeeze(mat2cell(p,size(p,1),size(p,2),repelem(1,size(p,3),1)))';
DT = delaunayTriangulation(p(:,1:2,1));
Faces = DT.ConnectivityList;
deformationStruct=triSurfaceDeformation_rewrited(Faces,pcell{1},pcell,1);
%%
figure; hold on; 
for iframe = 1:Nframe
plot(iframe,prctile(deformationStruct.Epc1{iframe},10,'all')*1e2,'.')
plot(iframe,prctile(deformationStruct.Epc1{iframe},90,'all')*1e2,'.')
plot(iframe,prctile(deformationStruct.Epc2{iframe},10,'all')*1e2,'.')
plot(iframe,prctile(deformationStruct.Epc2{iframe},90,'all')*1e2,'.')
end
%%
deformationStruct.FaceCentroids = cell(Nframe,1);
for ii = 1:Nframe
% compute face centroids
for iface=1:size(Faces,1)
deformationStruct.FaceCentroids{ii}(iface,:)=mean(pcell{ii}(Faces(iface,:),:));
end
end
%%
iframe = 30; 
figure; h_ax = gca; 
scatter3(deformationStruct.FaceCentroids{ii}(:,1),...
    deformationStruct.FaceCentroids{ii}(:,2),...
    deformationStruct.FaceCentroids{ii}(:,3),5,deformationStruct.Epc1{iframe}*1e2)
mean(deformationStruct.Epc1{iframe}*1e2)
set(h_ax(ax_idx),'XLim',[-15,15],'YLim',[-10,15],'ZLim',[-2,8],...
    'FontSize',12,'clipping','on','DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
colormap(h_ax(ax_idx),'jet');
colorbar(h_ax(ax_idx)); 
set(h_ax(ax_idx),'clim',[-10,10]);

%    fig = newfig('test interpolant 1');
%    Nrow = 1; Ncol = 1;
%    h_ax = subplot_ax(Nrow,Ncol); hold(h_ax,'on');
%    ax_idx = 1;
%    %     plot3(h_ax(ax_idx),PCcurrent(:,1),PCcurrent(:,2),PCcurrent(:,3),'r.')
%    mesh(xq,yq,vq1)
%    set(h_ax(ax_idx),'XLim',[-15,15],'YLim',[-10,15],'ZLim',[-2,8],...
%        'FontSize',12,'clipping','on','DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
%    colormap(h_ax(ax_idx),'jet');
%    set(h_ax(ax_idx),'clim',[-2,8]);


%    plot3(h_ax(ax_idx),xyJAI(:,1),xyJAI(:,2),vq1(k),'r.')
%    set(h_ax(ax_idx),...
%        'FontSize',12,'clipping','on','DataAspectRatio',[1 1 1],'Color','None');
%    colormap(h_ax(ax_idx),'jet');
%    set(h_ax(ax_idx),'clim',[-2,8]);