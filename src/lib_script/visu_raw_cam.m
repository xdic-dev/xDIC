
%% Visu settings 
subject = "S09"; 
material_nbr = 2;
trial = "007";
phase = "loading"; 
framejump = 1; 

[temp1,temp2] = import_vid(baseDataPath,...
    'subject',subject,...
    'material',material,...
    'trial',trial,...
    'stereopair',1,...
    'phase','all',...
    'framejump',framejump); 
[temp4,temp3] = import_vid(baseDataPath,...
    'subject',subject,...
    'material',material,...
    'trial',trial,...
    'stereopair',2,...
    'phase','all',...
    'framejump',framejump); 

%% Visu player 
fps_visu = actual_FS_vid/2; 
record_logic = false; 
player([temp1,temp2;temp4,temp3],'fps',fps_visu,'movie',record_logic);

% %%
% figure; 
% imhist(temp1(:,:,1)); 
% %%
% player(satur([temp2],'level',100),'movie',1,'fps',50); 
% %%
% player([temp4,temp3,temp2,temp1]); 
% 
% %%
% filename = "DIC3DPPresults_2Pairs_cum_v1"; 
% % multidicdata = load(fullfile('baseResultPath','analysis',subject,material,trial,phase,filename)); 
% ncorr1 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,"ncorr1")); 
% ncorr12 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,"ncorr12")); 
% ncorr2 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,"ncorr2")); 
% ncorr3 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,"ncorr3")); 
% ncorr43 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,"ncorr43")); 
% ncorr4 = load(fullfile(baseResultPath,'analysis',subject,material,trial,phase,"ncorr4")); 
% %%
% dy = 1216-1-200; 
% dx = 1936-850; 
% crop1 = [168,168+dx,175,175+dy]
% crop2 = [200,200+dx,150,150+dy]
% crop3 = [650,650+dx,150,150+dy]
% crop4 = [850,850+dx,1,1+dy]
% %%
% player(satur(temp4,'level',100))
% 
% 
% %%
% player(satur([temp1(crop1(3):crop1(4),crop1(1):crop1(2),:),...
%     temp2(crop2(3):crop2(4),crop2(1):crop2(2),:);...
%     temp4(crop4(3):crop4(4),crop4(1):crop4(2),:),...
%     temp3(crop3(3):crop3(4),crop3(1):crop3(2),:)...
%     ],'level',100),'movie',1,'fps',truefps/framejump); 
% %%
% roi1 = ncorr1.current_save(1).roi.mask; 
% roi2 = ncorr2.current_save(1).roi.mask; 
% roi3 = ncorr3.current_save(1).roi.mask; 
% roi4 = ncorr4.current_save(1).roi.mask; 
% %%
% 
% %%
% h = ncorr;
% h.set_ref(ncorr43.current_save(1).gs); 
% %%
% mask_out_3 = h.reference.roi.mask; 
% save(fullfile(baseResultPath,'fig','maskout_3.mat'),'mask_out_3'); 
% %%
% padding = 10; 
% %
% im3 = ncorr43.current_save(1).gs*255; 
% %im3(mask_out_3) = 0; 
% mask3 = ncorr3.reference_save.roi.mask; 
% im4 = ncorr43.current_save(2).gs*255; 
% %im4(mask_out_4) = 0; 
% mask4 = ncorr4.reference_save.roi.mask; 
% im1 = ncorr12.current_save(2).gs*255; 
% %im1(mask_out_1) = 0; 
% mask1 = ncorr1.reference_save.roi.mask; 
% im2 = ncorr12.current_save(1).gs*255; 
% %im2(mask_out_2) = 0; 
% mask2 = ncorr2.reference_save.roi.mask; 
% %
% Lymax = 0; 
% %
% [idx_y,idx_x] = find(ncorr3.reference_save.roi.mask == 1); 
% regionofinterest_y = min(idx_y)-3*padding:max(idx_y)+padding; Ly = length(regionofinterest_y); 
% regionofinterest_x = min(idx_x)-padding:max(idx_x)+4*padding; Lx = length(regionofinterest_x); 
% if Lymax < Ly
%     Lymax = Ly; 
% end
% %
% [idx_y,idx_x] = find(ncorr4.reference_save.roi.mask == 1); 
% regionofinterest_y = min(idx_y)-3*padding:max(idx_y)+padding; Ly = length(regionofinterest_y); 
% regionofinterest_x = min(idx_x)-padding:max(idx_x)+4*padding; Lx = length(regionofinterest_x); 
% if Lymax < Ly
%     Lymax = Ly; 
% end
% % 
% [idx_y,idx_x] = find(ncorr1.reference_save.roi.mask == 1); 
% regionofinterest_y = min(idx_y)-3*padding:max(idx_y)+padding; Ly = length(regionofinterest_y); 
% regionofinterest_x = min(idx_x)-4*padding:max(idx_x)+padding; Lx = length(regionofinterest_x); 
% if Lymax < Ly
%     Lymax = Ly; 
% end
% %
% [idx_y,idx_x] = find(ncorr2.reference_save.roi.mask == 1); 
% regionofinterest_y = min(idx_y)-3*padding:max(idx_y)+padding; Ly = length(regionofinterest_y); 
% regionofinterest_x = min(idx_x)-4*padding:max(idx_x)+padding; Lx = length(regionofinterest_x);
% if Lymax < Ly
%     Lymax = Ly; 
% end
% %
% [idx_y,idx_x] = find(ncorr3.reference_save.roi.mask == 1); 
% centery = floor((min(idx_y)-2*padding+max(idx_y)+padding)/2);
% regionofinterest_y = round(centery-Lymax/2):round(centery+Lymax/2);
% regionofinterest_x = min(idx_x)-padding:max(idx_x)+4*padding; 
% im3 = im3(regionofinterest_y,regionofinterest_x,:); 
% %mask3 = mask3(regionofinterest_y,regionofinterest_x); 
% %
% [idx_y,idx_x] = find(ncorr4.reference_save.roi.mask == 1); 
% centery = floor((min(idx_y)-2*padding+max(idx_y)+padding)/2);
% regionofinterest_y = round(centery-Lymax/2):round(centery+Lymax/2);
% regionofinterest_x = min(idx_x)-padding:max(idx_x)+4*padding; 
% im4 = im4(regionofinterest_y,regionofinterest_x,:); 
% %mask4 = mask4(regionofinterest_y,regionofinterest_x);
% % 
% [idx_y,idx_x] = find(ncorr1.reference_save.roi.mask == 1); 
% centery = floor((min(idx_y)-2*padding+max(idx_y)+padding)/2);
% regionofinterest_y = round(centery-Lymax/2):round(centery+Lymax/2);
% regionofinterest_x = min(idx_x)-4*padding:max(idx_x)+padding; 
% im1 = im1(regionofinterest_y,regionofinterest_x,:); 
% %mask1 = mask1(regionofinterest_y,regionofinterest_x); 
% %
% [idx_y,idx_x] = find(ncorr2.reference_save.roi.mask == 1); 
% centery = floor((min(idx_y)-2*padding+max(idx_y)+padding)/2);
% regionofinterest_y = round(centery-Lymax/2):round(centery+Lymax/2);
% regionofinterest_x = min(idx_x)-4*padding:max(idx_x)+padding; 
% im2 = im2(regionofinterest_y,regionofinterest_x,:); 
% %mask2 = mask2(regionofinterest_y,regionofinterest_x); 
% %
% 
% % figure; imshow(satur(im_highlight(im3,mask3,'factor',4),'level',100),[]);
% border = 20; 
% figure; 
% % imshow(satur([...
% %     im_highlight(im1,mask1,'factor',2),...
% %     zeros(Lymax+1,border)+255,...
% %     im_highlight(im2,mask2,'factor',2),...
% %     zeros(Lymax+1,border)+255,...
% %     im_highlight(im3,mask3,'factor',2),...
% %     zeros(Lymax+1,border)+255,...
% %     im_highlight(im4,mask4,'factor',2)...
% %     ],'level',100),[]); hold(gca,'on'); 
% 
% imshow(satur([im1,...
%     zeros(Lymax+1,border)+255,...
%     im2,...
%     zeros(Lymax+1,border)+255,...
%     im3,...
%     zeros(Lymax+1,border)+255,...
%     im4...
%     ],'level',100),[]); hold(gca,'on'); 
% ax = gca; 
% set(ax,'Color','None');
% 
% %%
% im1 = temp1(:,:,1:45*2); 
% im2 = temp2(:,:,1:45*2); 
% im3 = temp3(:,:,1:45*2); 
% im4 = temp4(:,:,1:45*2); 
% ylimit = 675;
% xlimit = 1075; 
% player(satur([im1((1:ylimit)+250,(1:xlimit)+0,:),im2((1:ylimit)+125,(1:xlimit)+230,:);...
%     im3((1:ylimit)+150,(1:xlimit)+850,:),im4((1:ylimit)+300,(1:xlimit)+750,:)],...
%     'level',100),...
%     'fps',25,...
%     'movie',1); 
% %%
% format = 'png';
% factor = 1; 
% %
% fig = figure; 
% set(fig,'Visible','off'); 
% imshow(satur([...
%     im_highlight(im1,mask1,'factor',factor)...
%     ],'level',100),[]); hold(gca,'on'); 
% ax = gca; 
% set(ax,'Color','None');
% set(gcf,'WindowState','maximized'); 
% exportgraphics(ax,fullfile(baseResultPath,'fig',sprintf('im1.%s',format)));hold(gca,'off'); 
% %
% fig = figure; 
% set(fig,'Visible','off'); 
% imshow(satur([...
%     im_highlight(im2,mask2,'factor',factor)...
%     ],'level',100),[]); hold(gca,'on'); 
% ax = gca; 
% set(ax,'Color','None');
% set(gcf,'WindowState','maximized'); 
% exportgraphics(ax,fullfile(baseResultPath,'fig',sprintf('im2.%s',format)));hold(gca,'off');
% %
% fig = figure; 
% set(fig,'Visible','off'); 
% imshow(satur([...
%     im_highlight(im3,mask3,'factor',factor)...
%     ],'level',100),[]); hold(gca,'on'); 
% ax = gca; 
% set(ax,'Color','None');
% set(gcf,'WindowState','maximized'); 
% exportgraphics(ax,fullfile(baseResultPath,'fig',sprintf('im3.%s',format)));hold(gca,'off');
% %
% fig = figure; 
% set(fig,'Visible','off'); 
% imshow(satur([...
%     im_highlight(im4,mask4,'factor',factor)...
%     ],'level',100),[]); hold(gca,'on'); 
% ax = gca; 
% set(ax,'Color','None');
% set(gcf,'WindowState','maximized'); 
% exportgraphics(ax,fullfile(baseResultPath,'fig',sprintf('im4.%s',format)));hold(gca,'off');
% %%
% % plot(idx_x(1),idx_y(1),'b.');
% % plot(idx_x(end),idx_y(end),'r.');
% % [~,idx] = max(idx_y); 
% % plot(idx_x(idx(end)),max(idx_y),'g.');
% % [~,idx] = min(idx_y); 
% % plot(idx_x(idx(end)),min(idx_y),'y.');
% %%
% max(max(uint8(ncorr43.current_save(1).gs)))





