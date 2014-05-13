% get coordinates from sector structure
sect1cord = sect1_2BH9.Coordinates
pdbid = sect1_2BH9.Pdb

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
        covMat(i,j) = 1/(numvect-1)*sum(diff_avg(:,i).*diff_avg(:,j));
    end
end
imagesc(covMat);
colorbar;

% calculate eigenvectors
[V,D] = eigs(covMat)


%% Plot data with eigenvectors

scatter3(sect1cord_x,sect1cord_y,sect1cord_z,20,'r','filled');
eigvect1 = V(:,1)*D(1,1);
eigvect2 = V(:,2)*D(2,2);
eigvect3 = V(:,3)*D(3,3);

sect_eigvectors = [eigvect1 eigvect2 eigvect3]

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
lambda_coordiantes = inv(V)*diff_avg';
figure
imagesc(lambda_coordiantes);
colorbar

%%
figure

scatter3(lambda_coordiantes(1,:), lambda_coordiantes(2, :),...
    lambda_coordiantes(3, :), 'filled')

title('Coordiantes of Sectors in Eigenspace', 'FontSize', 16)
xlabel('lambda 1', 'FontSize', 14) 
ylabel('lambda 2', 'FontSize', 14) 
zlabel('lamnda 3', 'FontSize', 14)


%%
% find the weight of each mode
sect_eigenvalues = diag(D);
variances = sect_eigenvalues/sum(sect_eigenvalues);
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
