function shaded_std_plot(x,y,varargin)
    %% Parse input arguments
    p = inputParser;
    p.addParameter('Color','k');
    p.addParameter('h_ax',gca);
    p.addParameter('FaceAlpha',0.5);
    p.addParameter('nstd',3);
    
    p.parse(varargin{:});
    c = p.Results.Color;
    h_ax = p.Results.h_ax;
    alpha = p.Results.FaceAlpha;
    nstd = p.Results.nstd;
    
    %
    mu = mean(y,2);
    dev = std(y,0,2);
    %
    ciplot(mu-nstd*dev,mu+nstd*dev,x,c,h_ax,alpha);
end