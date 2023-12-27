
%%% This code REQUIRES:
% ** The Optimization Toolboox
% intlinprog
% % using the following options
% % optimoptions(@intlinprog,'BranchRule','maxfun','Heuristics','none','Display','off');

% ** the seperate functions included
% find_leaky_set
% closure_comp


%% Computation for Z_2 (Q_3), takes approximately 10 seconds
disp('Optimal 2-leaky_set for Q_3:');
tic;
find_leaky_set(make_cube(3),2)
disp(['Computed in ', num2str(toc), ' seconds.',newline])

%% Computation for Z_3 (Q_4), takes approximately 200 seconds
disp('Optimal 3-leaky_set for Q_4');
tic;
Q=make_cube(4);
must_have = [];

% an open neighborhood of a vertex, excluding one neighbor, is a 3-leaky fort
for i=1:16
    z=zeros(16,1);
    z(i)=1;
    neighborhood=find(Q*z);
    must_have = [must_have; nchoosek(neighborhood,3)];
end
find_leaky_set(Q,3,must_have)
disp(['Computed in ', num2str(toc), ' seconds.',newline])

%% Verification for Z_3 (Q_5), takes approximately 10 seconds
disp('We know Q_4 is an optimal 0-leaky_set for Q_5; I will check Q_4 within Q_5 is an optimal 3-leaky forcing set.');

tic;
A = make_cube(5); %make Q_5
B=1; %trigger to end the verification in failure; 
% if a configuration of leaks is found where Q_4 does not work,
% B will set to 0;
ell = 3; %number of leaks 
w_index = 1:16; %indices for Q_4 within Q_5
F=[]; %list of potential forts; starts empty, forts will be appended upon failure

M = nchoosek(1:size(A,2),ell); %M is a set of possible leaks
    for i = 1:size(M,1) %go through each possible leak
        m = M(i,:);
        [b, v] = closure_comp(A, w_index, m); %find the closure complement
        % if the closure is not all vertices, b will be 0;
        if b == 0  %this loop converts v to a 0,1 vector for the F
            B=0;
            U = zeros(1,size(A,1)); %make a row of 0-1 indicator vector for the fort found.
            U(1,v) = 1;
            F = [F;U]; %append it to F
        end
    end

if B==1
    disp('Checked all configureations of 3 leaks; Q_4 (16 vertices) within Q_5 is a 3-leaky set!');
else
    disp('List of discovered forts:')
    F
end

disp(['Completed in ', num2str(toc), ' seconds.',newline])

%% Computation Z_1 for Grid Graphs, takes approximately 1200 seconds

for i=3:6
    disp('Computing 1-leaky forcing set for ',num2str(i),' x ',num2str(i), ' grid', newline);
    tic
    A=toeplitz([0,1,zeros(1,i-3),1],[0,1,zeros(1,i-3),1]); %creates a cycle
    AA=boxproduct(A,A);
    find_leaky_set(A,1)
    disp(['Completed ', num2str(i),' x ',num2str(i), ' grid in', num2str(toc), ' seconds.', newline])
end

%% Computation Z_1 for all csv in a folder, 
% takes approximately 10000 seconds for CubicGraphs

% requires you to select a folder

myDir = uigetdir; %gets directory
csvFiles = dir(fullfile(myDir,'*.csv'));
for csvFileName = csvFiles'
 A=readmatrix(csvFileName.name);
 disp(['Computing 1-leaky forcing set for ', csvFileName.name, ' ' , newline]);
 tic
 A=csvread(csvFileName.name);  %readthefile 

 if size(A,1)>2 && size(A,2) == 2 % if its an edgelist, make a matrix
     edgelist=A;
     A= zeros(max(max(edgelist)));
     for i=1:length(edgelist)
         A(edgelist(i,1),edgelist(i,2))=1;
         A(edgelist(i,2),edgelist(i,1))=1;
     end
 end

 if sum(sum(A))>0 %if A is not an empty graph
    find_leaky_set(A,1)
 end

    disp(['Completed ', csvFileName.name, ' in ', num2str(toc), ' seconds.', newline])
end

%% Computation Script for Z_1 for RandomGraphs, takes about 3000 seconds

z=0*(10:20); z=z'; 
t=0*(10:20); t=t';

numberoftrials=100;

for n=10:11
    currentTotals=[];
    currentTimes=[];

    for i=1:numberoftrials
        tic; u=length(find_leaky_set(randomGraph(n,1/4),1)); v=toc;
        currentTimes=[currentTimes v];
        currentTotals=[currentTotals u];
    end
    z(n)=mean(currentTotals);
    t(n)=mean(currentTimes);
end

numberoftrials=10;

for n=12:20
    currentTotals=[];
    currentTimes=[];

    for i=1:numberoftrials
        tic; u=length(find_leaky_set(randomGraph(n,1/4),1)); v=toc;
        currentTimes=[currentTimes v];
        currentTotals=[currentTotals u];
    end
    z(n)=mean(currentTotals);
    t(n)=mean(currentTimes);
end

%display the results
[z ; t]

%% auxillary functions to construct graphs

function [ G ] = make_cube( d, graphoutq )
%Makes graph G with dimension d; G will be outputed as a matrix unless 
%graphoutq is defined and nonzero

if nargin<2
    graphoutq=0;
end

%d = max(1,floor(d)) %makes d make sense
A = make_path(2);
B = make_path(2);
Z = boxproduct(A,B);

if d==1
    A=G;
    return
end

for n = 3:d
    A = Z;
    Z = boxproduct(A,B);
end

if graphoutq
    G = graph(Z);
else
    G=Z;
end   
end

function [ A ] = boxproduct( G, H )
%Creates the adjacency matrix of the box product of G and H
A = kron(G,eye(length(H))) + kron(eye(length(G)),H);
end

function [A] = randomGraph(n,p)
% Makes an Erdos-Reyni random graph, output as an adjacency matrix
A=triu(rand(n,n)<p); A=A+A';
end

function [ A ] = make_path( n )
%Creates the adjacency matrix of a path on n vertices

A = zeros(n);
for i = 2:n-1
    A(i,i-1)=1;
    A(i,i+1)=1;
end

A(1,2)=1;
A(n,n-1)=1;
end

