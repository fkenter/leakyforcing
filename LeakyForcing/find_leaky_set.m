function [w_opt, F] = find_leaky_set(A, ell, must_have)
% finds a leaky forcing set of graph with adjacency matrix A, ell leaks,
%   containing necessarily the must_have vertices.

% input: A is the adjacency matrix of the graph
%   ell is the number of leaks
%   must_have is a list of the vertices (index) that must be in the leaky
%       forcing set
% output: an optimal t leaky forcing set as a list of indicies

options=optimoptions(@intlinprog,'BranchRule','maxfun','Heuristics','none','Display','off');
%options for the 

% default values
if nargin < 3
    must_have = [];
end

if nargin == 1
    ell = 0;
end

B = 0; %B is a breaker variable; will turn 1 when loop is to be broken.
F = ones(1,size(A,1)); % F is the set of forts.
S = sum(A,1); %S is the degrees

%Adds to the list of forts- vertices with degree \leq number of leaks
for i = 1:size(A,2)
    U = zeros(1,size(A,1)); %temp array to build a fort of uncolored vertices
    if S(i) <= ell
        U(1,i) = 1;
        F = [F;U]; %add the uncolored vertices to the list of forts.
    end

end

%Adds to the list of forts- vertices in must_have must be in the zero forcing set
for i = 1:size(must_have,1)
    U = zeros(1,size(A,1));
    U(must_have(i,:)) = 1;
    F = [F;U];
end

fortcounter=1000;

while B == 0 %while there is not an optimal solution found,
    B=1; % B will turn to 0 if there is a configuration of new leaks that won't work;
    f = ones(size(A,2),1); %all-one vector for objective function
    r = size(F,1); %size of F for the constraint
    %counter to print progress for larger problems

    %run the integer program. ** Requires optimization package ** 
    w = intlinprog(f, 1:size(f,1), -F, -ones(r,1),[],[], zeros(size(f,1),1),...
        ones(size(f,1),1),options);

    % finds the index for each non zero entry of the optimal solution
    w_index = find(w);
    w_index = w_index';

    M = nchoosek(1:size(A,2),ell); %M is a set of possible leaks
    for i = 1:size(M,1) %go through each possible leak
        m = M(i,:);
        [b, v] = closure_comp(A, w_index, m); %find the closure complement
        % b is an indicator of whether or not the closure covered all
        % vertices.

        if b == 0  %this loop converts v to a 0,1 vector for the F
            B=0;
            U = zeros(1,size(A,1)); %make a row of 0-1 indicator vector for the fort found.
            U(1,v) = 1;
            F = [F;U]; %append it to F
            if size(F,1) > fortcounter
                disp([num2str(fortcounter), ' forts found so far']);
                fortcounter = fortcounter + 1000;
            end
        end
    end
end
w_opt = w_index; %return the last optimal set of colored vertices.
end