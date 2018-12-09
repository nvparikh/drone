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

global p_degree time_per_segments total_element derivative_cache num_segments;
persistent way_points parameters times derivatives pDegree;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 2 %First time call
    p_degree = 5;
    pDegree = p_degree;
    total_element = p_degree + 1;
    derivative_cache = cell(p_degree,1);
    for index=1:p_degree
        derivative_cache{index} = getPositionDerivatives(index-1);
    end
    derivatives = derivative_cache; %persisting
    %calculate euler distance between waypoints and adjust time accordingly
    d = waypoints(:,2:end) - waypoints(:,1:end-1);
    d0 = 2 * sqrt(d(1,:).^2 + d(2,:).^2 + d(3,:).^2);
    num_segments = size(waypoints,2) -1;
    time_for_1_distance = 1;
    time_per_segments = cumsum(d0*time_for_1_distance)
    times = time_per_segments; %to retail times for each segments
    
    parameters = zeros(num_segments*total_element, size(waypoints,2));
    for i=1:size(waypoints,1) %iterate for xyz
        [A_EQ_segment_boundary, B_eq_segment_boundary] = getContinuityConstraintForSegmentBoundary();
        [A_EQ_position, B_eq_position] = getAbsoluteConstraints(waypoints(i,:)', 1, [0 time_per_segments]);
        remaining_constraints = num_segments*total_element - (size(A_EQ_segment_boundary,1)+size(A_EQ_position,1));
        num_seg_1_extra_const = ceil(remaining_constraints/2);
        num_last_seg_extra_const = floor(remaining_constraints/2);
        [A_EQ_remaining_seg_1, B_eq_remaining_seg_1] = getPerSegmentConstraints(zeros(num_seg_1_extra_const,1), 2:2+num_seg_1_extra_const-1, 1, 0);
        [A_EQ_remaining_seg_last, B_eq_remaining_seg_last] = getPerSegmentConstraints(zeros(num_last_seg_extra_const,1), 2:2+num_last_seg_extra_const-1, num_segments, time_per_segments(num_segments));
        A_EQ=[A_EQ_segment_boundary;A_EQ_position;A_EQ_remaining_seg_1;A_EQ_remaining_seg_last];
        B_eq=[B_eq_segment_boundary;B_eq_position;B_eq_remaining_seg_1;B_eq_remaining_seg_last];
        parameters(:,i) = A_EQ\B_eq;
    end
    parameters
    way_points = waypoints;
else
    if(t > times(end))
        t = times(end);
    end
    seg_no = find(times >= t,1);
    pos = derivatives{1}(t);
    vel = derivatives{2}(t);
    acc = derivatives{3}(t);
    waypoint_polynomial_index_start = (seg_no-1)*(pDegree+1)+1;
    waypoint_polynomial_indexes = waypoint_polynomial_index_start:waypoint_polynomial_index_start+pDegree;
    desired_state.pos = [pos*parameters(waypoint_polynomial_indexes,1);pos*parameters(waypoint_polynomial_indexes,2);pos*parameters(waypoint_polynomial_indexes,3)];
    desired_state.vel = [vel*parameters(waypoint_polynomial_indexes,1);vel*parameters(waypoint_polynomial_indexes,2);vel*parameters(waypoint_polynomial_indexes,3)];
    desired_state.acc = [acc*parameters(waypoint_polynomial_indexes,1);acc*parameters(waypoint_polynomial_indexes,2);acc*parameters(waypoint_polynomial_indexes,3)];
    desired_state.yaw = 0;
    desired_state.yawdot = 0;
    % desired_state.pos = zeros(3,1);
% desired_state.vel = zeros(3,1);
% desired_state.acc = zeros(3,1);
% desired_state.yaw = 0;
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
    v = zeros(total_length, 1);
    v(start_pos:start_pos+length(x)-1) = x;
end

function [A_EQ, B_eq] = getContinuityConstraintForSegmentBoundary()
    global p_degree time_per_segments derivative_cache total_element num_segments;

    matlabOffset = 1;
    maxSmoothDerivativeAtBoundary = p_degree - 1;

    A_EQ = zeros((num_segments - 1)*maxSmoothDerivativeAtBoundary, total_element*num_segments);
    B_eq = zeros((num_segments - 1)*maxSmoothDerivativeAtBoundary, 1);

    for index= 1:num_segments-1
        time = time_per_segments(index);
        start = (index -1)*total_element + matlabOffset;
        for smoothDerivative =1:maxSmoothDerivativeAtBoundary
            constraintRow = derivative_cache{smoothDerivative}(time);

            %The way we make sure derivative is same at segment boundary is by
            %[polynomial_derivative, - polnomial_derivate,...] *
            %[params1,params2,...]' = [0...]. Thats why B_eq is all zero
            %Also we generate each row for as many smooth derivative like
            %position, velocity, acceleration etc. we want
            A_EQ((index - 1)*maxSmoothDerivativeAtBoundary + smoothDerivative, start:start+total_element*2-1) = [constraintRow, -constraintRow];
        end
    end
end

function [A_EQ, B_eq] = getPerSegmentConstraints(values, derivatives, segment_no, time)
    global derivative_cache total_element num_segments;
    size_of_values = size(values,1);
    A_EQ = zeros(size_of_values, num_segments*total_element);
    B_eq = values;
    start_index = (segment_no-1)*total_element+1;
    for i = 1:size(values,1)
        A_EQ(i, start_index:start_index+total_element-1) = derivative_cache{derivatives(i)}(time);
    end
end

function [A_EQ, B_eq] = getAbsoluteConstraints(values, nth_derivative, times)
    global derivative_cache total_element;
    size_of_values = size(values,1);
    A_EQ = zeros(size_of_values, (size_of_values-1)*total_element);
    B_eq = values;
    for i = 1:size_of_values
        start_index_segment = getSegStartFromConstraintChain(i);
        A_EQ(i, start_index_segment:start_index_segment+total_element-1) = derivative_cache{nth_derivative}(times(i));
    end
end

function start_index_segment = getSegStartFromConstraintChain(i)
    global total_element
    seg_no = (i-1);
    if i==1
       seg_no = 1;
    end
    start_index_segment = (seg_no-1)*total_element+1;
end




