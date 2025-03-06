function varargout = fit_xyplane_pc(points,varargin)
    %% Parse input arguments
    p = inputParser;
    p.addRequired('points');
    p.addParameter('margin',0);

    p.parse(points,varargin{:});
    margin = p.Results.margin;
    
    %% Computation  
    % Design Matrix
    DM = [points(:,1),points(:,2),ones(length(points(:,1)),1)]; 
    % Estimate Parameters
    B = DM\points(:,3);
%     
%     figure; hold(gca,'on');
    % Plane function (x,y) : return z coordinate  
    fun_planefit = @(x,y) B(2)*y + (B(3)+margin);%B(1)*x + 
    % Point below fitted plane 
    idx_in = fun_planefit(points(:,1),points(:,2)) > points(:,3);
%     scatter3(points(:,1),points(:,2),points(:,3),5,[1 1 1].*0.5,'filled'); 
%     scatter3(points(idx_in,1),points(idx_in,2),points(idx_in,3),10,'g','filled');
% %     scatter3(points(:,1),points(:,2),fun_planefit(points(:,1),points(:,2)),5,[1 1 1].*0.5,'filled'); 
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
    
%     B = robustfit(DM,points(:,3),[],[],true)
%     figure; hold(gca,'on');
%     % Plane function (x,y) : return z coordinate  
%     fun_planefit = @(x,y) B(3)*y + (B(1)+0.1);%B(1)*x + 
%     % Point below fitted plane 
%     idx_in = fun_planefit(points(:,1),points(:,2)) > points(:,3);
%     scatter3(points(:,1),points(:,2),points(:,3),5,[1 1 1].*0.5,'filled'); 
%     scatter3(points(idx_in,1),points(idx_in,2),points(idx_in,3),10,'g','filled');
% %     scatter3(points(:,1),points(:,2),fun_planefit(points(:,1),points(:,2)),5,'k','filled'); 
%     set(gca,'DataAspectRatio',[1 1 1],'Color','None','View', [0 -90]);
    
    %% outputs arguments
    varargout{1} = idx_in; 
    varargout{2} = fun_planefit; 
    
end