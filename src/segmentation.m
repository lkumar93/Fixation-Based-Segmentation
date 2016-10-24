%
% THIS IS AN IMPLEMENTATION OF THE PAPER ACTIVE SEGMENTATION WITH
% FIXATION BY AJAY MISHRA ET AL AT ICCV 2009
%
% COPYRIGHT BELONGS TO THE AUTHOR OF THIS CODE 
%
% AUTHOR : LAKSHMAN KUMAR
% AFFILIATION : UNIVERSITY OF MARYLAND, MARYLAND ROBOTICS CENTER
% EMAIL : LKUMAR93@UMD.EDU
% LINKEDIN : WWW.LINKEDIN.COM/IN/LAKSHMANKUMAR1993
%
% THE WORK (AS DEFINED BELOW) IS PROVIDED UNDER THE TERMS OF THE MIT LICENSE
% THE WORK IS PROTECTED BY COPYRIGHT AND/OR OTHER APPLICABLE LAW. ANY USE OF
% THE WORK OTHER THAN AS AUTHORIZED UNDER THIS LICENSE OR COPYRIGHT LAW IS PROHIBITED.
% 
% BY EXERCISING ANY RIGHTS TO THE WORK PROVIDED HERE, YOU ACCEPT AND AGREE TO
% BE BOUND BY THE TERMS OF THIS LICENSE. THE LICENSOR GRANTS YOU THE RIGHTS
% CONTAINED HERE IN CONSIDERATION OF YOUR ACCEPTANCE OF SUCH TERMS AND
% CONDITIONS.
%


close all
clear all

addpath(genpath('../maxflow-v3.0/'))

%% Read image
imName = 'Images/DSC01527';
fix_x = 209;   % fixation point
fix_y = 145;   % fixation point

im = imread([imName, '.png']);
im2 = zeros(size(im));
im = im2double(im);
edgeGrad = load([imName, '_grad.mat']);  %% Use precomputed Berkeley edges
edgeGrad = edgeGrad.edgeGrad;

%   im,             input image (must be double)
%   edgeGrad,       image probabilistic edge map
%   fix_x,          fixation point x coordinate
%   fix_y,          fixation point y coordinate

%% Energy function parameters

% Binary weights
nju = 5;                    % Binary weight exponent constant
k = 20;                     % Binary weight for two pixels with zero edge probability
lambda = 1000;              % Importance of binary weights

% Unary weights
D  = 1e100;                 % Fixation point unary weight

%unaryColourWeight = 1;      % Importance of colourmodel unary weights

foreground = 2;             % Object
background = 1;             % Background

%% Image parameters
[yRes, xRes] = size(edgeGrad);
nPix = xRes*yRes;

%% Compute binary weights

E = edges4connected(yRes,xRes);                                             % Indices of adjacent pixels (Potts model)

% The algorithm proceeds as follows
% 1. Average edge probability at adjacent edges
% 2. Edge where at least one of the pixels belongs to the edge map 
% 3. Edges where none of the pixels belong to the edge map (assign them to k)
% 4. Calculate the distance of each edge from the fixation point
% 5. Weights are the inverse of the distance from the fixation point
% 6. Normalize the weights to have maximum of 1
% 7. Construct unary weights for image boundary and fixation point 
% 8. Perform min-cut using max flow function
% 


%Get Size of Edges Connected

[m1,n1] = size(E);

AvgEdgeGrad = zeros(size(E,1),1);

%Get average edge probability at adjacent edges by adding elements in
%the 2 columns and dividing the sum by 2
for i = 1: m1

    AvgEdgeGrad(E(i,1)) = (edgeGrad(E(i,1)) + edgeGrad(E(i,2)))/2;    
    
end

%Initialize Binary Weights
BinaryWeights = zeros(size(E,1),1);

%Assign BinaryWeight as k if Edge Probability for the pixel is zero or to
%exp(-nju*(EdgeProbability(p)+EdgeProbability(q))/2)
for i = 1: m1
    
    if   AvgEdgeGrad(E(i,1)) ~= 0

        BinaryWeights(i) = exp(-nju*(AvgEdgeGrad(i)))  ;
        
    else
        BinaryWeights(i) = k; 
    end
        
end

%Initialize the distances for the pixels from fixation point
Distances = zeros(size(E,1),1);

%Compute the midpoint of neighbouring pixels and find the distances between
%the midpoint and the fixation point

for i = 1: m1
    
    [y1,x1] = ind2sub([yRes xRes],E(i,1));
    [y2,x2] = ind2sub([yRes xRes],E(i,2));
    xmid = (x1 + x2)/2 ;
    ymid = (y1 + y2)/2;
    Distances(i) =   sqrt((ymid-fix_y)^2 + (xmid-fix_x)^2);
end

%Initialize the weights to multiplied with the already computed binary
%weights
Weights = zeros(size(E,1),1);

%Weights for corresponding pixels are given by inverse of its distance from
%the fixation point
for i =1 :m1
    Weights(i) = 1/(Distances(i));
end

%Find the maximum among all the weights
MaxWeight = max(Weights);

%Normalize the weights to be between 0 and 1
for i =1 :m1
    Weights(i) = Weights(i)/MaxWeight;
end

%Compute the modified binary weights by multiplying it with the weights
for i =1 :m1
    BinaryWeights(i) = BinaryWeights(i)*Weights(i);
end

%Compute the spare matrix corresponding to Binary Weights

A = sparse(E(:,1),E(:,2),BinaryWeights,nPix,nPix,4*nPix);


%% Compute Unary Weights


%Initialize the Unary Weights
UnaryWeights = zeros(prod(size(edgeGrad)),2);

%Assign unary weight D to borders of the image and the fixation point

%Assign unary weight D to Left most Column of the image
for i =1 : yRes
    UnaryWeights(i,background) = D;
end

%Assign unary weight D to Right most Column of the image
for i = (prod(size(edgeGrad))-yRes +1): -1 : prod(size(edgeGrad))
    UnaryWeights(i,background) = D;
end

%Assign unary weight D to Top Most Row of the image

for i = 0:(xRes-1)
    UnaryWeights(i*yRes+1,background) = D;
end

%Assign unary weight D to Bottom Most Row of the image

for i = 1:xRes
    UnaryWeights(i*yRes,background) = D;
end

%Assign unary weight D to the fixation point
ind = sub2ind([yRes,xRes],fix_y,fix_x);
UnaryWeights(ind,foreground) =  D;

%Create a sparse matrix of UnaryWeights

T = sparse(UnaryWeights);


%% Perform Min-Cut and Display Segmented Image

%Perform graph min-cut using  maxflow function
[flow,labels] = maxflow(A,T);

%Reshape the labels to fit the image size
labels = reshape(labels, [yRes,xRes]);

%Display the segmented object

%If the label for corresponding pixel is 1 , display that pixel

for i=1 : yRes
    for j=1 :xRes
        for k = 1 :3
        
            if(labels(i,j) == 1)

                im2(i,j,k) =im(i,j,k);

            end
            
        end
    end
end


figure;

image(im2);
title('Segmented Image');


%% Explanation Of The Paper

% EXPLANATION
% Before starting the actual segmentation process, we should generate an 
% edge map of the given input images wherein the actual boundary edges of 
% the object are bright, while the internal edges that are due to textures  
% and intensity are dim. In order to do this first , the images are passed 
% through  the Berkeley edge detector which will give the initial edge 
% map, in which edges due textures are removed. But still it will have 
% some internal non-boundary edges. To suppress these strong internal 
% edge segments and to strengthen the boundary edges, we break the 
% initial edge map into straight line segments and select rectangular
% regions of certain width at a certain distance on both sides of the 
% edges. Since we have stereo images, we can find the average disparity 
% inside these rectangles. The difference in the average disparity is the
% measure of likelihood of the segment to be the depth boundary. 
% Let change in disparity be denoted as ?d. 
% Then the brightness of the edge pixels on the segment is changed.
%  
% After modifying the edge as per the equation above, the image borders 
% are also added as edges so as to ensure that objects that are partially 
% visible in the images are enclosed. However the intensity of these edges
% are kept low such that they are not preferred over real edges.
% 
% 
% Now that we have the probabilistic edge map, which is in Cartesian 
% co-ordinates. We need to convert this to polar co-ordinates. For this we 
% need a 2d Gaussian Kernel function which is applied to all the edge pixels.

%  
%  
% 
% Now the polar edge map is calculated by sampling the W(x,y) whose 
% intensities are scaled to be between 0 and 1. And each row of edge map
% corresponds to different angles(0 to 360), while each column of the edge 
% map corresponds to different radii (0 to r_max). This will have a graph
% like structure , where each pixel ,apart from those in the first and 
% last columns of the polar edge map, are connected to four immediate
% neighbours. Here the first and last rows of the polar edge map are also
% connected to each other. All these pixels have to be labeled either as
% 0 if the pixel is inside the edge or 1 if the pixel is outside the 
% edge. Pixels in the first column are labeled as 0, as it contains the 
% fixation point itself and the pixels in the last column are labeled as 1,
% as they surely would be outside the edges. In order to label the remaining
% pixels in the other column , we use the graph cut algorithm on the energy 
% function  , so as to find the labeling function which corresponds to the 
% minimum energy.
% 
%  
%  
% 
% 
% Where the costs for assigning the labels are given as below
% U_p (l_p=1)=D and U_p (l_p=0)=0  
%  
% The resulting binary segmentation splits the polar edge map into two
% parts: left part which contains all the pixels inside the edge and the 
% right part which contains all the pixels outside the edge. However this 
% segmentation does not take into account the color of the pixels . 
% The color information at  a pixel in the polar map is obtained by 
% interpolating the RGB value at the corresponding sub-pixel location
% in the Cartesian Space	.
% 
% The color distributions of the pixels inside and outside can be
% considered as normalized  3d histograms with 10 bins along each channel
% (Red,Green & Blue). Now the data/cost term for all the pixels except those
% in the first and last columns becomes 
%  
% We again use graph cut algorithm to minimize the modified energy function
% with the new data term. The segmentation results improve after introducing
% the color information in the energy formulation.
% 







    
