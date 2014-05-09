function [sector_pca ] = sector_pca(sector_structure)
% will perform pca on defined protein sectors
%   Output from Sector.m will give centroids of a coevolving sector.  Here
%   we will do PCA on ONE of the sectors. The output of sector_pca is a 3D
%   scatterplot of the protein sector and the eigenvectors, deconvolution
%   heat map, and bar graph of variences by mode number.


% get coordinates from sector structure
sect1cord = sector_structure.Coordinates
pdbid = sector_structure.Pdb

% extract coodinates
sect1cord_x = sect1cord(1, :)
sect1cord_y = sect1cord(2, :)
sect1cord_z = sect1cord(3, :)

% define matrix for PCA
pca_mat = zeros(length(sect1cord), 3);
pca_mat(:, 1) = sect1cord_x;
pca_mat(:, 2) = sect1cord_y;
pca_mat(:, 3) = sect1cord_z;


% calculate mean

avg = mean(pca_mat)

% differance avg

numvect = size(pca_mat, 1); %number of vectors
diff_avg = pca_mat-repmat(avg,numvect,1);
scatter(diff_avg(:,1), diff_avg(:,2),5, 'filled')

% covariance matrix
ndim = size(pca_mat,2); 
for i=1:ndim 
    for j=1:ndim
        cov(i,j) = 1/(numvect-1)*sum(diff_avg(:,i).*diff_avg(:,j));
    end
end
imagesc(cov);
colorbar;

% calculate eigenvectors
[V,D] = eigs(cov)
%% Plot data with eigenvectors

scatter3(sect1cord_x,sect1cord_y,sect1cord_z,20,'r','filled');
eigvect1 = V(:,1)*D(1,1);
eigvect2 = V(:,2)*D(2,2);
eigvect3 = V(:,3)*D(3,3);

hold on;
plot3([avg(1)-eigvect1(1) avg(1)+eigvect1(1)], [avg(2)-eigvect1(2) avg(2)+eigvect1(2)],...
    [avg(3)-eigvect1(3) avg(3)+eigvect1(3)], 'Color',  [0 0 1], 'LineWidth', 3);

plot3([avg(1)-eigvect2(1) avg(1)+eigvect2(1)], [avg(2)-eigvect2(2) avg(2)+eigvect2(2)],...
    [avg(3)-eigvect2(3) avg(3)+eigvect2(3)], 'Color',  [0 0 .7], 'LineWidth', 3);

plot3([avg(1)-eigvect3(1) avg(1)+eigvect3(1)], [avg(2)-eigvect3(2) avg(2)+eigvect3(2)],...
    [avg(3)-eigvect3(3) avg(3)+eigvect3(3)], 'Color',  [0 0 .3], 'LineWidth', 3);
hold off
axis equal

legend('centroids', 'eigenvector 1', 'eigenvector 2', 'eigenvector 3')
xlabel('x coordinates', 'FontSize', 14)
ylabel('y coordinates', 'FontSize', 14)
zlabel('z coodinates', 'FontSize', 14)
title('PCA of eigenvectors', 'FontSize', 16)
%%
%decompose the vectors
amounts = inv(V)*diff_avg';
figure
imagesc(amounts);
colorbar

% find the weight of each mode
eigenvalues = diag(D);
variances = eigenvalues/sum(eigenvalues);
for i=1:3
    fprintf('The fraction of the variance from mode %d is %g\n',i,variances(i));
end
figure
bar(variances)
xlabel('mode number', 'FontSize', 14)
ylabel('distance (angstroms)', 'FontSize', 14)
title('Fraction of variance of sector by mode number', 'FontSize', 16)
% 
end

