function [G] = graph_of_2d_array(unit);
% creates an unweighted graph of a 2d array with starting node and ending nodes, connected to
% the first and last rows of the array, respectively (for purposes of Dijkstra's algorithm)
 
% create adjacency matrix with nodes = elems of matrix, and connect all
% numerical neighbours and up/down and diagonally (per heatmap matrix)
elems = numel(unit);
dims = size(unit);
row_scale = dims(1)/dims(2)*100;
col_scale = dims(2)/dims(1);
graph_conn = zeros(elems) + diag([repmat(row_scale,elems-1,1)],1) + di-ag([repmat(row_scale,elems-1,1)],-1) + diag([repmat(row_scale*col_scale,elems-50,1)],50) + diag([repmat(row_scale*col_scale,elems-50,1)],-50) + di-ag([repmat(row_scale*4,elems-49,1)],49) + diag([repmat(row_scale*4,elems-49,1)],-49) + diag([repmat(row_scale*4,elems-51,1)],51) + di-ag([repmat(row_scale*4,elems-51,1)],-51);
 
% remove edges joining opposite sides of matrix ('cut the cylinder')
unit_dims = size(unit);
for K = unit_dims(1):unit_dims(1):elems
	graph_conn(K,K+1) = 0; % remove end-to-end joining, e.g. if unit_dims(1)= 50, then remove 50-51
    graph_conn(K+1,K) = 0;
    try
        graph_conn(K-49,K-50) = 0; % (mirror)
        graph_conn(K-50,K-49) = 0;
		graph_conn(K+1,K-50) = 0; % remove end-to-end diagonals, 
		% e.g. if unit_dims(1)= 50, then remove 101-50
        graph_conn(K-50,K+1) = 0;
        graph_conn(K,K-49) = 0; % (mirror)
        graph_conn(K-49,K) = 0;
    catch
	    fprintf('Error cutting ‘cylinder’ graph for this unit')
    end
end
 
graph_conn(1,50) = 0; graph_conn(50,1) = 0;
graph_conn = graph_conn(1:elems,1:elems); % remove phantom nodes
 
% create start and end nodes
start_node = elems+1;
end_node = elems+2;
graph_conn(start_node,1:unit_dims(1)) = 1; % start node
graph_conn(1:unit_dims(1),start_node) = 1;
graph_conn(end_node,elems-unit_dims(1):elems) = 1; % end node
graph_conn(elems-unit_dims(1):elems,end_node) = 1;
 
% make directed graph
G = digraph(graph_conn);
 
end

sets = {'L4-Hartmann-sham','L4-Basic-sham','L4-Basic-tbi','L2-Hartmann-sham','L2-Hartmann-tbi','L2-Basic-sham','L2-Basic-tbi'};
 
for set = 1:length(sets)
    data = sets{1,set};
    dashes = strfind(data,'-');
    layer = data(1:2);
    health = data(dashes(2)+1:end);
    type = data(dashes(1)+1:dashes(2)-1);
    
    events = eval([layer,health,'EventSummary_',type]);
    name_col = set;
    figDir = ['figures\',layer,'\'];
 
    for K = 1:length(events(:,1,1))
        unit = events(K,:,:);
        unit = squeeze(unit);
        unit = unit*20;
        len = numel(unit);
        
        try
            % create graph and weight edges
            [G] = graph_of_2d_array(unit);
            plot(G) % optional
            weights = abs(unit/max(max(unit)) - 1) * 100; % normalise
            weights(weights==0) = 1; % turn 0s into 1s
            weights = reshape(weights,1,len);
 
            % assign weights to graph
            dest_nodes_list = G.Edges.EndNodes(:,2);
            for node = 1:len
                G.Edges.Weight(dest_nodes_list==node) = weights(node);
            end
 
            % use Dijkstra's algorithm to find shortest path
            [path1,d1] = shortestpath(G,len+1,len+2,'method','positive');
 
            % find 'second' shortest path
			weights(path1(2:end-1)) = path1(2:end-1) + 1000; % add arb large weights to nodes from 'first' shortest path
            for node = 1:len
                G.Edges.Weight(dest_nodes_list==node) = weights(node);
            end
			[path2,d2] = shortestpath(G,len+1,len+2,'method','positive'); % find path
 
            % find 'third' shortest path
			weights(path2(2:end-1)) = path2(2:end-1) + 1000; % add arb large weights to nodes from 'second' shortest path
            for node = 1:len
                G.Edges.Weight(dest_nodes_list==node) = weights(node);
            end
			[path3,d3] = shortestpath(G,len+1,len+2,'method','positive'); % find path
 

            paths1{set,K,:} = path1;
            paths2{set,K,:} = path2;
            paths3{set,K,:} = path3;
 
            dists1(set,K,:) = d1;
            dists2(set,K,:) = d2;
            dists3(set,K,:) = d3;
        catch
            fprintf(['Error for: set ',set,', unit ',num2str(K)])
        end
        
        
    end
    
    clear weights
    
end
 
dijkstra_analysis = {{paths1,paths2,paths3},{dists1,dists2,dists3}};