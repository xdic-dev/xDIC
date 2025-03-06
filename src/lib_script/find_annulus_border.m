function idx_IN_annulus = find_annulus_border(PCnow,CAnow,varargin)

%% Parse input arguments
p = inputParser;
p.addRequired('PCnow'); 
p.addRequired('CAnow'); 
p.addParameter('borderWidth',1);

p.parse(PCnow,CAnow,varargin{:});
border_width_annulus = p.Results.borderWidth;

%% find point inside and annulus around the contact area 
idx_IN_annulus = false(size(PCnow,1),1);
if ~isnan(CAnow)
    [~,centerEl,Rx,Ry,theta_radians] = myfitellipse(CAnow(:,1),CAnow(:,2)); 
    ratio = Ry/Rx;
    Rxlarge = Rx+border_width_annulus; 
    Rxsmall = Rx-border_width_annulus; 
    Rylarge = Rxlarge*ratio; 
    Rysmall = Rxsmall*ratio; 

    %for small ellipse at the contact creation 
    if Rxsmall<0 || Rysmall<0
        Rxsmall = 0; 
        Rysmall = 0; 
    end

    %fit ellipse around 
    CAfitlarge = create_ellipse(centerEl(1),centerEl(2),Rxlarge,Rylarge,theta_radians);
    CAfitsmall = create_ellipse(centerEl(1),centerEl(2),Rxsmall,Rysmall,theta_radians);    

    %find idx inside the border
    idx_IN_curv = inpolygon(PCnow(:,1),PCnow(:,2),CAfitlarge(:,1),CAfitlarge(:,2));
    idx_IN_small = inpolygon(PCnow(:,1),PCnow(:,2),CAfitsmall(:,1),CAfitsmall(:,2));
    idx_IN_annulus = idx_IN_curv&~idx_IN_small;

%     figure; hold on;
%     scatter3(PCnow(:,1),PCnow(:,2),PCnow(:,3),1,[1 1 1].*0.5);
%     scatter3(PCnow(idx_IN_annulus,1),PCnow(idx_IN_annulus,2),PCnow(idx_IN_annulus,3),5,'r');
%     plot3(CAnow(:,1),CAnow(:,2),CAnow(:,3),'b','LineWidth',2);
%     
%     xlim(gca,[-15,15]);
%     ylim(gca,[-10,25]);
%     zlim(gca,[-2,15]);
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [90 0]);
%     hold off;

end
end