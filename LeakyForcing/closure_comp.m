function [b, v] = closure_comp(A, color, leaks)

%Determines whether or not 'color' (a list of vertices) is a leaky forcing set on a graph A
% that has leaks on vertices in list 'leaks' (a list of vertices)
% Output: b is a 0-1 boole that indicates whether or not it colors all
% vertices
%   v is a list of uncolored vertices that result, e.g., a fort.

if nargin < 3 %if leaks are not specified, there are none.
    leaks = [];
end

C = zeros(size(A)); %matrix for storage of blue neighbors as a 0-1 array
% e.g. C(i,j) is 1 is j is blue.
U = zeros(size(A)); %matrix for storage of white neighbors as a 0-1 array
% e.g. U(i,j) is 1 is j is white.

%Creating initial C and U matrices
for i = 1:size(A,2)
    for j = 1:size(A,2)
        if A(i,j) == 1 && ismember(j, color) == 1
            C(i,j) = 1;
        end

        if A(i,j) == 1 && ismember(j,color) == 0 
            U(i,j) = 1;
        end
    end
end

%checks to see if a colored vertex has only one uncolored neighbor
for k = 1:size(A,2) 
    %The number of iterations of the rule are at most
    %the number of vertices
    new_color = [];
    for i = 1:size(A,2)
        if sum(U(i,:)) == 1 && ismember(i,color) == 1 && ismember(i,leaks) == 0
            z = find(U(i,:));
            new_color = [new_color z]; %stores which colored vertices have one neighbor uncolored
        end
        
    end
    
    color = [color new_color];
    color = unique(color);
    
    %changes C and U matrices
    for i = 1:size(color,2)
        for j = 1:size(A,2)
            if A(j,color(i)) == 1
                U(j,color(i)) = 0;
                C(j,color(i)) = 1;
            end
        end
    end
end


% sets return values
if size(color,2) == size(A,2) % if completely colored, return b=1 and the empty set
    b = 1;
    v = [];
else  % otherwise, return b=0 and the set of uncolored vertices.
    b = 0;
    v = 1:size(A,2);
    v(color) = []; 
end
end