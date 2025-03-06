function C = compute_covariance_matrix2(points)
tic
logic_valid = ~isnan(points(:,1));
points = points(logic_valid,:); 
mean_xyz = mean(points,1); 
C = (mean_xyz - points)'*(mean_xyz-points)/size(points,1); 
toc
end