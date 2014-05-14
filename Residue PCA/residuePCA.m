function [eigVect, eigVal, aaEigenBase] = residuePCA ( aaDBcount, normalize, graph )
	% PCA on normalized data
	if normalize == 1
		% Transpose to have aa as column and
		% sectors as row
		aaCount = aaDBcount';

		% Mean occurence for each aa
		avg = mean(aaCount,1);

		% Degree of freedom (nDim) and number of sectors
		[numSector, nDim] = size(aaCount);

		% Covariance matrix
		diff_avg = aaCount-repmat(avg,numSector,1);
		for i=1:nDim
	    	for j=1:nDim
	        	cov(i,j) = 1/(numSector-1)*sum(diff_avg(:,i).*diff_avg(:,j));
	       	end
	    end



	else
	    % PCA on raw data
	    aaCount = aaDBcount';

		% Degree of freedom (nDim) and number of sectors
		[numSector, nDim] = size(aaCount);

		% Covariance matrix
		for i=1:nDim
	    	for j=1:nDim
	        	cov(i,j) = 1/(numSector-1)*sum(aaCount(:,i).*aaCount(:,j));
	       	end
	    end
	end

    % Eigen values and vectors
    % The eigen values are along the diagonal of eigVal
    % They are not ordered, I think.
    [V,D] = eig(cov);
    eigVect = V;
    eigVal = D;

    % Coordinates in the eigen base:
    aaEigenBase = inv(V)*aaCount';

 	eigenvalues = diag(D);
	variances = eigenvalues/sum(eigenvalues);
	variances = floor(1000*variances)/10;
	for i=1:20
    	fprintf('Mode %d explains %g percents of the variance.\n',i,variances(i));
	end



	if graph == 1
	    % Projection in 2D
	    aaProjection2dim = aaEigenBase(end-1:end,:);
	    figure
	    scatter(aaProjection2dim(2,:),aaProjection2dim(1,:));
	    title('Projection on the 2 highest modes');
	    xlabel('PCA 20');
	    ylabel('PCA 19');

	    % Projection in 3D
	    aaProjection3dim = aaEigenBase(end-2:end,:);
	    figure
	    scatter3(aaProjection3dim(3,:),aaProjection3dim(2,:),aaProjection3dim(1,:));
	    title('Projection on the 3 highest modes');
	    xlabel('PCA 20');
	    ylabel('PCA 19');
	    zlabel('PCA 18');
	end

end