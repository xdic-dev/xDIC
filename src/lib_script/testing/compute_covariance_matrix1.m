function C = compute_covariance_matrix1(points)
tic;
logic_valid = ~isnan(points(:,1));
points = points(logic_valid,:); 
mean_xyz = mean(points,1); 
C = zeros(3,3); 
for ii = 1:3
    for jj = 1:3
        for kk = 1:size(points,1)
            C(ii,jj) = C(ii,jj) + (mean_xyz(ii) - points(kk,ii))*(mean_xyz(jj)-points(kk,jj));
        end
    end
end
C = C./size(points,1); 
toc;
end