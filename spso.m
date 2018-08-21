% Ethan Fahnestock 2018
% Soumya D. Mohanty, Jul 2018: modified interface to match GENERICPSO/PSO
% 
% Argument  descriptions 
% nDim: number of dimensions in the search space
% 
% ff_pointer: matlab handle to the function
% 
% psoParams: input parameters
%     .pop_size: number of  particles: Should be set to 40
%     .maxSteps: number of iterations of the main loop
%     .c: value of the c constant in the G calculation: Should be set to
%     0.5+ln(2)
%     .w: inertia of velocity: Should be set to 1/(2ln(2))
% 
% returns bestPos, (best position in range 0-1)
%         bestFit, best fitness found
function varargout = spso(ff_pointer, nDim, psoParams, outputParams)

bestFit = Inf; % tracks the best fitness value
bestPos = zeros(1,nDim); %tracks the best position

k_nbrs = 3; %% number of neighbors excluding self
neighborhoods = randi([1,psoParams.pop_size], psoParams.pop_size, k_nbrs);

%% initialize particles
pop(psoParams.pop_size).pos = rand(1,nDim); %initialze last member to preallocate
for i = 1:psoParams.pop_size
    pop(i).pos = rand(1,nDim);
    pop(i).vel = rand(1,nDim) -  pop(i).pos;
    pop(i).cur_fit = [];
    pop(i).lbest_fit = Inf;
    pop(i).pbest_fit = Inf;
    pop(i).pbest = pop(i).pos;
    pop(i).lbest = []; % just initialze vector
end

% run main loop
for i = 1:psoParams.maxSteps
    % update fitness values
    improved = false;  % tracks if improvements were made in a loop iteration 
    for  j = 1:psoParams.pop_size
       pop(j).cur_fit = ff_pointer(pop(j).pos);
       if pop(j).cur_fit < pop(j).pbest_fit % update pbest
            pop(j).pbest = pop(j).pos;
            pop(j).pbest_fit = pop(j).cur_fit;
       end
       
       if pop(j).cur_fit < bestFit % if we have found a new gbest
           bestFit = pop(j).cur_fit;
           bestPos = pop(j).pos;
           improved = true;
       end
       
    end % end particle fitness loop
    
    if ~improved % if we did not improve after a run randomize the neighborhoods again
        neighborhoods = randi([1,psoParams.pop_size], psoParams.pop_size, k_nbrs);
    end
    
    %lbest update loop
    for j = 1:psoParams.pop_size
        bst_nbhr_fit = Inf;
        bst_nbhr_index = []; %initalize scalar w/ no value
        for k = 1:k_nbrs % for every neighbor
            nbr_index = neighborhoods(j, k);
            if pop(nbr_index).cur_fit < bst_nbhr_fit % set to best found neighbor
                bst_nbhr_fit = pop(nbr_index).cur_fit;
                bst_nbhr_index = nbr_index;
            end
        end  % end  neighbor loop
        if pop(j).cur_fit < bst_nbhr_fit %check particle itself, as it is always a neighbor of itself
            bst_nbhr_fit = pop(j).cur_fit;
            bst_nbhr_index = j;
        end
            
         if bst_nbhr_fit < pop(j).lbest_fit % if we have found a new lbest fitness value
             pop(j).lbest_fit = bst_nbhr_fit;
             pop(j).lbest = pop(bst_nbhr_index).pos;
         end
        
    end % end lbest update loop
    
    % velocity, position, and constriction loop
    for j = 1:psoParams.pop_size
        %calculate g
        if pop(j).pbest == pop(j).lbest 
            G = pop(j).pos + psoParams.c * (pop(j).pbest - pop(j).pos)/2.0;
        else % if lbest != pbest
            G = pop(j).pos + psoParams.c * (pop(j).pbest + pop(j).lbest - 2.0 * pop(j).pos) / 3.0;
        end %end lbest == pbest
        
        %create  x_prime,  random point within hypersphere
        x_prime = normrnd(0,1,[1,nDim]);
        x_prime = x_prime / norm(x_prime) * rand() * norm(G-pop(j).pos); %  scale random hypersphere point to random radius
        
        x_prime = x_prime + G; % move center to G
        
        %update velocity
        pop(j).vel = psoParams.w * pop(j).vel  + x_prime - pop(j).pos;
        % update position
        pop(j).pos = pop(j).pos + pop(j).vel;
        
%         %apply confinement
%         for k = 1:nDim
%             if abs(pop(j).vel(k)) > 0.5
%                 pop(j).vel(k) = 0.5 * sign(pop(j).vel(k));
%             end
%         end % end confinement
    end % end velocity, pos, update loops
end % end main loop

%Optional Output
if nargout 
    returnData = struct('totalSteps',[],...
        'totalFuncEvals',[],...
        'bestLocation',bestPos,...
        'bestSNR',bestFit);
    varargout{1}=returnData;
    if nargout > 1
        varargout{2}=[];
    end
end
end % end function