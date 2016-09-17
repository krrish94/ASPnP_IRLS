function bestInd = pick_best_solution(U, u, R_est, t_est)
% PICK_BEST_SOLUTION  Used when one of the cases of ASPnP returns multiple
% solutions, i.e., (R,t) pairs. This returns the index of the solution with
% the least reprojection error.
%
% The symbols U, u, R_est, t_est are consistent with the nomenclature
% adopted in the ASPnP_IRLS code.


% By default, assume that index 1 is the best solution
bestInd = 1;

% Initialize the least reprojection error to inf
leastReprojError = inf;

% For each solution to be examined
for i = 1:size(t_est,2)
    % Compute the reprojection error
    proj = R_est(:,:,i)*U + t_est(:,i)*ones(1,size(u,2));
    proj = proj ./ repmat(proj(3,:),3,1);
    reprojErr = abs(proj(1:2,:) - u);
    reprojErr = sum(sum(reprojErr));
    % Update if it is better than the currently selected solution
    if reprojErr < leastReprojError
        bestInd = i;
        leastReprojError = reprojErr;
    end
end


end