function fillgray(ax,x1,x2)
    ylim = get(ax,'YLim'); 
    fill(ax,[x1 x2 x2 x1], [ylim(1) ylim(1) ylim(2) ylim(2)],[1 1 1].*0.0,'FaceAlpha',0.3)
end