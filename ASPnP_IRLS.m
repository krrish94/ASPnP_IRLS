function [R0, t0, w_final] = ASPnP_IRLS(U0, u0, K, w_prior)
% ASPNP_IRLS  Implements the ASPnP method from the paper "ASPnP: An 
% Accurate and Scalable Solution to the Perspective-n-Point Problem" in 
% the IEICE Trans. on Information and Systems.
% The paper isn't publicly available (it's on ResearchGate, though). But,
% almost the entire algorithm can be found in the paper "Revisiting the PnP
% Problem: A Fast, General and Optimal Solution".
%
% I make a small modification to incorporate an IRLS loop for outlier
% rejection.
%
% Inputs:
%   U0: 3D points (3-by-n matrix)
%   u0: 2D projections of the corresponding 3D points (2-by-n matrix)
%   K: intrinsic camera matrix
%   w_prior: a prior on the weights (based on the viewpoint)


% Assume normalized image coordinates if K is not passed as input. Else,
% perform normalization (pre-multiply by inv(K)).
if nargin < 3
    u = u0(1:2,:);
else
    u = K\[u0(1:2,:);ones(1,size(u0,2))];
    u = u(1:2,:);
end

% Number of points
n = size(U0, 2);

% If the prior weights are not specified, assume them to be uniform
if nargin < 4
    w_prior = ones(n,1);
end

% Normalize 3D points
% [U C3D] = normalize_3D_points(U0);

% Inhomogeneous coordinates (if U0 has 1's in the last row, i.e., if it is
% a 4-by-n matrix, drop the last row)
U = U0(1:3,:);


% (3D) Matrix to hold the rotation estimate(s) obtained from various cases
C_est = [];
% Matrix to hold the translation estimate(s) obtained from various cases
t_est = [];

% We take all possible solutions for the rotation matrix R (using the
% proposed non-unit quaternion parametrization) and split them into 4
% 'independent' cases.


%% Testing IRLS formulation with case 1

% Normalize the prior weights
w_prior = exp(w_prior) ./ sum(exp(w_prior));

% Number of iterations for which IRLS is to be run
numIters = 5;

% Initialize IRLS with the uniform weight solution

% Run case 1
[R_est_1, t_est_1] = ASPnP_IRLS_Case1_V2(U, u, w_prior);

% Hack for now. Replace this with code that retains the best (R,t). (Issue closed now. Hack replaced.)
% R_est_1 = R_est_1(1:3,1:3,1);
% t_est_1 = t_est_1(:,1);

% The above function frequently returns more than one solution. In such a
% situation, iterate over all solutions and pick the best one.
bestInd = pick_best_solution(U, u, R_est_1, t_est_1);
R_est_1 = R_est_1(:,:,bestInd);
t_est_1 = t_est_1(:,bestInd);

% % The following lines of code needn't be used. Their functionality is 
% % now included in the pick_best_solution function.
% if size(t_est_1,2) > 0
%     bestInd = 1;
%     leastReprojErr = inf;
%     % Initialize the solution with the first pair of (R,t)
%     soln_R_est_1 = R_est_1(:,:,1);
%     soln_t_est_1 = t_est_1(:,1);
%     % For every solution pair
%     for i = 1:size(t_est_1,2)
%         % Compute the reprojection error to each keypoint
%         proj = R_est_1(:,:,i)*U + t_est_1(:,i)*ones(1,size(u,2));
%         proj = proj ./ repmat(proj(3,:),3,1);
%         reprojErrs = abs(proj(1:2,:) - u);
%         reprojErrs = sum(sum(reprojErrs)');
%         % Update if it is better than the currently selected solution
%         if reprojErrs < leastReprojErr
%             bestInd = i;
%             leastReprojErr = reprojErrs;
%         end
%     end
%     % Update R_est_1 and t_est_1
%     R_est_1 = R_est_1(:,:,bestInd);
%     t_est_1 = t_est_1(:,bestInd);
% end


% Run IRLS iterations
for i = 1:numIters
    
    % Compute the reprojection error to each keypoint
    proj = R_est_1*U + t_est_1*ones(1, size(u,2));
    proj = proj ./ repmat(proj(3,:),3,1);
    reprojErrs = abs(proj(1:2,:) - u);
    reprojErrs = sum(reprojErrs)';
    
    % Re-weight each keypoint, ensuring that the prior weights also take
    % part in the re-weighting process
    prior_mult = 0.3;
    w = 1./(reprojErrs + 0.1);
    w = exp(w) ./ sum(exp(w));
    % Alternative weight updates
    w = prior_mult*w_prior + (1-prior_mult)*w;
    w = exp(w) ./ sum(exp(w));
    
    % Re-run case 1, with the new weights
    [R_est_1, t_est_1] = ASPnP_IRLS_Case1_V2(U, u, w);
    % If more than one solution is returned, pick a solution with the least
    % reprojection error
    bestInd = pick_best_solution(U, u, R_est_1, t_est_1);
    R_est_1 = R_est_1(:,:,bestInd);
    t_est_1 = t_est_1(:,bestInd);
    
end

% Store final weights (to be returned)
w_final = w;

% Final (R,t)
R0 = R_est_1;
t0 = t_est_1;


%% Original method proposed in the paper

% % Case 1: When a ~= 0, without loss of generality, we can rescale a, b, c,
% % d implicitly, such that a = 1. Then, there are only 3 remaining
% % parameters b, c, d in R. This is same as the Cayley parametrization used
% % in DLS (Direct Least Squares) for PnP.
% [R_est_1, t_est_1] = ASPnP_IRLS_Case1_V2(U, u);
% % If a rotation-translation pair is estimated, save it
% if size(t_est_1,2) > 0
%     C_est = cat(3,C_est,R_est_1);
%     t_est = [t_est t_est_1];
% end
% % Case 2: When a = 0, b ~= 0, by letting b = 1, it is possible to reduce
% % the rotation matrix parametrization.
% [R_est_2, t_est_2] = ASPnP_IRLS_Case2_V2(U, u);
% % If a rotation-translation pair is estimated, save it
% if size(t_est_2,2) > 0
%     C_est = cat(3,C_est,R_est_2);
%     t_est = [t_est t_est_2];
% end
% % Case 3: When a = 0, b = 0, c ~= 0, we can let c = 1
% [R_est_3, t_est_3] = ASPnP_IRLS_Case3_V2(U, u);
% % If a rotation-translation pair is estimated, save it
% if size(t_est_3,2) > 0
%     C_est = cat(3,C_est,R_est_3);
%     t_est = [t_est t_est_3];
% end
% % Case 4: When a = 0, b = 0, c = 0, d ~= 0, we can let d = 1
% [R_est_4, t_est_4] = ASPnP_IRLS_Case4_V2(U, u);
% % If a rotation-translation pair is estimated, save it
% if size(t_est,4) > 0
%     C_est = cat(3,C_est,R_est_4); t_est = [t_est t_est_4];
% end


%% Select best solution

% % Among all these solutions, choose the one with the smallest reprojection
% % error. (What do you do if you have multiple solutions with the same
% % reprojection error???)
% 
% % Create a vector of ones, equal to the number of solutions (usually 4)
% index = ones(1,size(t_est,2));
% % For each solution
% for i = 1:size(t_est,2)
%     % Project the 3D points to the image, using the ith solution
%     proj = C_est(:,:,i)*U+t_est(:,i)*ones(1,size(u,2));
%     % If points are not in front of the camera, ignore the index (set it to
%     % zero so that it is automatically ignored in the later code)
%     if min(proj(3,:)) < 0
%         index(i) = 0;
%     end
%     % Compute reprojection error for the current solution
%     proj = proj./repmat(proj(3,:),3,1);
%     err = proj(1:2,:)-u;
%     error(i) = sum(sum(err.*err));
% end
% % Choose the valid solution with the smallest reprojection error. Valid
% % solution implies that all points must be in front of the camera.
% error0 = min(error(index>0));
% 
% % If no solution is valid, return an infinite cost
% if isempty(error0)
%     R0 = [];
%     t0 = [];
%     cost = inf;
%     return;
% end


%% Set output variables

% % Using 5% as tolerance to detect multiple solutions
% rr = find((error < error0*1.05).*(index>0));
% R0 = C_est(:,:,rr);
% t0 = t_est(:,rr);
% % Cost (reprojection error) of the solution
% cost = error(rr);


end
