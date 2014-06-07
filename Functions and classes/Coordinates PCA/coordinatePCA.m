function [eigValues, eigVectors, coordinatesEigenBase] = coordinatePCA(sector, showplots)
% will perform pca on defined protein sectors
%   Output from Sector.m will give centroids of a coevolving sector.  Here
%   we will do PCA on ONE of the sectors. The output of sector_pca is a 3D
%   scatterplot of the protein sector and the eigenvectors, deconvolution
%   heat map, and bar graph of variences by mode number.


    % get coordinates from sector structure
    pca_mat = sector.Coordinates';

    % calculate mean
    avg = mean(pca_mat);

    % differance avg
    numvect = size(pca_mat, 1); %number of vectors;
    diff_avg = pca_mat-repmat(avg,numvect,1);

    % covariance matrix
    ndim = size(pca_mat,2); 
    for i=1:ndim ;
        for j=1:ndim;
            covMat(i,j) = 1/(numvect-1)*sum(diff_avg(:,i).*diff_avg(:,j));
        end;
    end;

    % calculate eigenvectors
    [eigVectors, eigValues] = eigs(covMat);
    eigValues = diag(eigValues);


    %decompose the vectors
    coordinatesEigenBase = inv(eigVectors)*diff_avg';



    %%
    % find the weight of each mode
    variances = eigValues/sum(eigValues);
    for i=1:3;
        fprintf('The fraction of the variance from mode %d is %g\n',i,variances(i));
    end;
    %%
    if showplots == 1
        figure;
        bar(variances);
        xlabel('mode number', 'FontSize', 14);
        ylabel('distance (angstroms)', 'FontSize', 14);
        title('Fraction of variance of sector by mode number', 'FontSize', 16);

        figure
        scatter3(diff_avg(:,1), diff_avg(:,2), diff_avg(:,3), 5, 'filled');
        title('Scatter plot of the residue in the x,y,z planes')

        figure
        imagesc(covMat);
        colorbar;
        title('Covariation matrix')

        %% Plot data with eigenvectors
        figure
        eigvect1 = eigVectors(:,1)*eigValues(1);
        eigvect2 = eigVectors(:,2)*eigValues(2);
        eigvect3 = eigVectors(:,3)*eigValues(3);
        scatter3(pca_mat(:,1) ,pca_mat(:,2) , pca_mat(:,3), 20,'r','filled');
        hold on;
        plot3([avg(1)-eigvect1(1) avg(1)+eigvect1(1)], [avg(2)-eigvect1(2) avg(2)+eigvect1(2)],...
        [avg(3)-eigvect1(3) avg(3)+eigvect1(3)], 'Color',  [0 0 1], 'LineWidth', 3);

        plot3([avg(1)-eigvect2(1) avg(1)+eigvect2(1)], [avg(2)-eigvect2(2) avg(2)+eigvect2(2)],...
        [avg(3)-eigvect2(3) avg(3)+eigvect2(3)], 'Color',  [0 0 .7], 'LineWidth', 3);

        plot3([avg(1)-eigvect3(1) avg(1)+eigvect3(1)], [avg(2)-eigvect3(2) avg(2)+eigvect3(2)],...
        [avg(3)-eigvect3(3) avg(3)+eigvect3(3)], 'Color',  [0 0 .3], 'LineWidth', 3);
        hold off;
        axis equal;

        legend('centroids', 'eigenvector 1', 'eigenvector 2', 'eigenvector 3');
        xlabel('x coordinates', 'FontSize', 14);
        ylabel('y coordinates', 'FontSize', 14);
        zlabel('z coodinates', 'FontSize', 14);
        title('PCA of eigenvectors', 'FontSize', 16);

        variances = eigValues/sum(eigValues);
        for i=1:3;
            fprintf('The fraction of the variance from mode %d is %g\n',i,variances(i));
        end;

        % Projection in eigenbase
        figure;
        scatter3(coordinatesEigenBase(1,:), coordinatesEigenBase(2, :),...
            coordinatesEigenBase(3, :), 'filled');

        title('Coordiantes of Sectors in Eigenspace', 'FontSize', 16);
        xlabel('eigenvector 1', 'FontSize', 14);
        ylabel('eigenvector 2', 'FontSize', 14);
        zlabel('eigenvector 3', 'FontSize', 14);

    end

end

