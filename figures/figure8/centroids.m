function cents = centroids(samples); 
% returns the centroid location 

num = size(samples,1);
addpath ~/Desktop/softwares/triangular_dislocation/

cents = zeros(num,3); 
parfor i = 1: num
    model = samples(i,:); 
    [~,~,X,Y,Z] = Planar_to_tri(model,0);
    
    cents(i,:) = [mean(X); mean(Y); mean(Z)];  
end

    