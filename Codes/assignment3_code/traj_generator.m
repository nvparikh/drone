function [ desired_state ] = traj_generator(t, state, waypoints)
% TRAJ_GENERATOR: Generate the trajectory passing through all
% positions listed in the waypoints list
%
% NOTE: This function would be called with variable number of input arguments.
% During initialization, it will be called with arguments
% trajectory_generator([], [], waypoints) and later, while testing, it will be
% called with only t and state as arguments, so your code should be able to
% handle that. This can be done by checking the number of arguments to the
% function using the "nargin" variable, check the MATLAB documentation for more
% information.
%
% t,state: time and current state (same variable as "state" in controller)
% that you may use for computing desired_state
%
% waypoints: The 3xP matrix listing all the points you much visited in order
% along the generated trajectory
%
% desired_state: Contains all the information that is passed to the
% controller for generating inputs for the quadrotor
%
% It is suggested to use "persistent" variables to store the waypoints during
% the initialization call of trajectory_generator.


%% Example code:
% Note that this is an example of naive trajectory generator that simply moves
% the quadrotor along a stright line between each pair of consecutive waypoints
% using a constant velocity of 0.5 m/s. Note that this is only a sample, and you
% should write your own trajectory generator for the submission.

global p_degree time_per_segments total_element derivative_cache num_segments num_coordinates parameters_length;
persistent way_points parameters times derivatives pDegree noCoordinates;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 2 %First time call
    p_degree = 6;
    num_coordinates = size(waypoints,1);
    noCoordinates = num_coordinates;
    pDegree = p_degree;
    total_element = p_degree + 1;
    derivative_cache = cell(p_degree,1);
    for index=1:p_degree
        derivative_cache{index} = getPositionDerivatives(index-1);
    end
    derivatives = derivative_cache; %persisting
    %calculate euler distance between waypoints and adjust time accordingly
    d = waypoints(:,2:end) - waypoints(:,1:end-1);
    d0 = sqrt(d(1,:).^2 + d(2,:).^2 + d(3,:).^2);
    num_segments = size(waypoints,2) -1;
    time_for_1_distance = 2;
    time_per_segments = cumsum(d0*time_for_1_distance)
    times = time_per_segments; %to retail times for each segments
    snapDerivativeIndex = 5;
    parameters_length = num_coordinates*num_segments*total_element;

    [A_EQ_segment_boundary, B_eq_segment_boundary] = getContinuityConstraintForSegmentBoundary(3)
    [A_EQ_position, B_eq_position] = getAbsoluteConstraints(waypoints', 1, [0 time_per_segments])
    %remaining_constraints = num_segments*total_element - (size(A_EQ_segment_boundary,1)+size(A_EQ_position,1));
    num_seg_1_extra_const = pDegree-1; %position, velocity, accleration = 0 at beggining
    num_last_seg_extra_const = pDegree-1; %position, velocity, accleration = 0 at end
    [A_EQ_remaining_seg_1, B_eq_remaining_seg_1] = getPerSegmentConstraints(zeros(num_seg_1_extra_const,num_coordinates), 2:2+num_seg_1_extra_const-1, 1, 0)
    [A_EQ_remaining_seg_last, B_eq_remaining_seg_last] = getPerSegmentConstraints(zeros(num_last_seg_extra_const,num_coordinates), 2:2+num_last_seg_extra_const-1, num_segments, time_per_segments(num_segments))
    A_EQ=[A_EQ_segment_boundary;A_EQ_position;A_EQ_remaining_seg_1;A_EQ_remaining_seg_last];
    B_eq=[B_eq_segment_boundary;B_eq_position;B_eq_remaining_seg_1;B_eq_remaining_seg_last];
    %parameters(:,i) = A_EQ\B_eq;
    %Solving constraint optimization problem to find parameters
    %objectiveFun can be passed, but in this code it is hard coded to
    %integration(square(snap)). Integration can be done analytically in this
    %case and used. But we are numerically integrating square(snap)
    objFunction = @(params) integral(@(t) objectiveFunValue(t, snapDerivativeIndex, params), 0,time_per_segments(end));

    %Col vector for all segements [param_seg1, param_seg2...]'
    initialParams = 10*rand(parameters_length, 1);
    options = optimoptions('fmincon','Algorithm','sqp');
    [params, fval, exitflag, output] = fmincon(objFunction, initialParams, [], [], A_EQ, B_eq,[],[],[],options);
    disp(output);
    parameters = params;

    way_points = waypoints;
else
    if(t > times(end))
        t = times(end);
    end
    seg_no = get_segment_no(t);
    seg_index_start = (seg_no-1)*(pDegree+1);

    desired_state.pos = zeros(noCoordinates, 1);
    desired_state.vel = zeros(noCoordinates, 1);
    desired_state.acc = zeros(noCoordinates, 1);

    for coor = 1:noCoordinates
        coor_start = get_coordinate_start(coor) + seg_index_start;
        param_indeces = coor_start:coor_start+pDegree;
        desired_state.pos(coor) = derivatives{1}(t)*parameters(param_indeces);
        desired_state.vel(coor) = derivatives{2}(t)*parameters(param_indeces);
        desired_state.acc(coor) = derivatives{3}(t)*parameters(param_indeces);
    end
    desired_state.yaw = 0;
    desired_state.yawdot = 0;
end
end

function dpdx = getPositionDerivatives(numDerivatives)
    global p_degree total_element;
    x = [];
    syms x;
    positionWOParams = sym('x', [1  total_element]);
    for degree=0:p_degree
        positionWOParams(degree + 1) = power(x, degree);
    end
    dpdx = matlabFunction(diff(positionWOParams, numDerivatives));
    %matlabFunction greatly improves speed
end
    
function v = pad_with_zero(x, start_pos, total_length)
    v = zeros(1,total_length);
    v(start_pos:start_pos+length(x)-1) = x;
end

function start_index = get_coordinate_start(coor_no)
    global total_element num_segments
    start_index = (coor_no-1)*num_segments*total_element + 1;
end

function [A_EQ, B_eq] = getContinuityConstraintForSegmentBoundary(maxSmoothness)
    global time_per_segments derivative_cache total_element num_segments num_coordinates parameters_length;

    A_EQ = zeros((num_segments - 1)*num_coordinates*maxSmoothness, parameters_length);
    B_eq = zeros((num_segments - 1)*num_coordinates*maxSmoothness, 1);

    iter = 1;
    for coor = 1:num_coordinates
        coor_start = get_coordinate_start(coor); 
        for index= 1:num_segments-1
            time = time_per_segments(index);
            start = (index -1)*total_element + coor_start;
            for smoothDerivative =1:maxSmoothness
                constraintRow = derivative_cache{smoothDerivative}(time);

                %The way we make sure derivative is same at segment boundary is by
                %[polynomial_derivative, - polnomial_derivate,...] *
                %[params1,params2,...]' = [0...]. Thats why B_eq is all zero
                %Also we generate each row for as many smooth derivative like
                %position, velocity, acceleration etc. we want
                A_EQ(iter, start:start+total_element*2-1) = [constraintRow, -constraintRow];
                iter = iter + 1;
            end
        end
    end
end

%values=  no_value*cordinate
function [A_EQ, B_eq] = getPerSegmentConstraints(values, derivatives, segment_no, time)
    global derivative_cache total_element num_coordinates parameters_length;
    length_of_values = size(values,1);
    A_EQ = zeros(length_of_values*num_coordinates, parameters_length);
    B_eq = reshape(values, [], 1);
    iter = 1;
    for coor = 1:num_coordinates
        start_index = (segment_no-1)*total_element+get_coordinate_start(coor);
        for i = 1:length_of_values
            A_EQ(iter, start_index:start_index+total_element-1) = derivative_cache{derivatives(i)}(time);
            iter = iter + 1;
        end
    end
end

%values=  no_value*cordinate
function [A_EQ, B_eq] = getAbsoluteConstraints(values_per_seg, nth_derivative, times)
    global derivative_cache total_element num_coordinates num_segments parameters_length;
    total_points_per_coor = num_segments + 1;
    A_EQ = zeros(total_points_per_coor, parameters_length);
    B_eq = reshape(values_per_seg, [], 1);
    iter = 1;
    for coor = 1:num_coordinates
        for i = 1:total_points_per_coor
            start_index_segment = getSegStartFromConstraintChain(i) + get_coordinate_start(coor);
            A_EQ(iter, start_index_segment:start_index_segment+total_element-1) = derivative_cache{nth_derivative}(times(i));
            iter = iter + 1;
        end
    end
end

function start_index_segment = getSegStartFromConstraintChain(waypointIndex)
    global total_element
    seg_no = (waypointIndex-1);
    if waypointIndex==1
       seg_no = 1;
    end
    start_index_segment = (seg_no-1)*total_element;
end

function segNo = get_segment_no(time)
    global time_per_segments
    segNo = find(time_per_segments >= time,1);
end

function objVals = objectiveFunValue(times, objectiveDerivative, parameters)
    global derivative_cache total_element num_coordinates parameters_length
    A_EQ = zeros(size(times,2)*num_coordinates, parameters_length);
    iter = 1;
    for index=1:size(times,2)
        t = times(index);
        polynomial_value = derivative_cache{objectiveDerivative}(t);
        segNo = get_segment_no(t);
        segStartIndex = (segNo-1)*total_element;
        for coor = 1:num_coordinates
            startIndex = segStartIndex + get_coordinate_start(coor);
            A_EQ(iter, startIndex:startIndex+total_element-1) =  polynomial_value;
            iter = iter + 1;
        end
    end
    objVals = sum(reshape((A_EQ*parameters).^2,num_coordinates, size(times,2)),1);
end





