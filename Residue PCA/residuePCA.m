function [eigVect, eigVal, aaEigenBase] = residuePCA ( aaDBcount, ...
												 normalize, graph )
%RESIDUEPCA calculates the covariance matrix of the matrix of aa
%occurences and returns its eigen vectors and eigenvalues, as well
%as the matrix of aa occurences in the eigen base.
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

	elseif normalize == 2
		% Transpose to have aa as column and
		% sectors as row
		aaCount = aaDBcount';

		% Frequency of amino acid from Uniprot
		aaFreq = [ 8.25 ; 5.53 ; 4.06 ; 5.45 ; 1.37 ; 3.93 ; ...
				6.75 ; 7.07 ; 2.27 ; 5.96 ; 9.66 ; 5.84 ; 2.42 ; ...
				3.86 ; 4.70 ; 6.56 ; 5.34 ; 1.08 ; 2.92 ; 6.87 ]'/100;

		% Degree of freedom (nDim) and number of sectors
		[numSector, nDim] = size(aaCount);

		% Normalizing by removing the average frequenc of
		% amino acids in uniprot
		normAaCount = aaCount - repmat(aaFreq, numSector, 1);

		% Covariance matrix
		for i=1:nDim
	    	for j=1:nDim
	        	cov(i,j) = 1/(numSector-1)*sum(normAaCount(:,i).*normAaCount(:,j));
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
	    title(['Projection of the matrix of amino acid' char(10) ...
	    		'counts on the 2 highest modes']);
	    xlabel('PCA 20');
	    ylabel('PCA 19');

	    % Projection in 3D
	    aaProjection3dim = aaEigenBase(end-2:end,:);
	    figure
	    scatter3(aaProjection3dim(3,:),aaProjection3dim(2,:),aaProjection3dim(1,:));
	    title(['Projection of the matrix of amino acid' char(10) ...
	    		'counts on the 3 highest modes']);
	    xlabel('PCA 20');
	    ylabel('PCA 19');
	    zlabel('PCA 18');
	end

end