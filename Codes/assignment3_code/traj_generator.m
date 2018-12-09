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
    p_degree = 6;
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
    time_for_1_distance = 1;
    time_per_segments = cumsum(d0*time_for_1_distance)
    times = time_per_segments; %to retail times for each segments
    corridorWidths = [2.43];
    noCorridorConstraints = 20;
    snapDerivativeIndex = 5;
    
    parameters = zeros(num_segments*total_element, size(waypoints,2));
    for i=1:size(waypoints,1) %iterate for xyz
        [A_EQ_segment_boundary, B_eq_segment_boundary] = getContinuityConstraintForSegmentBoundary(3);
        [A_EQ_position, B_eq_position] = getAbsoluteConstraints(waypoints(i,:)', 1, [0 time_per_segments]);
        %remaining_constraints = num_segments*total_element - (size(A_EQ_segment_boundary,1)+size(A_EQ_position,1));
        num_seg_1_extra_const = 3; %position, velocity, accleration = 0 at beggining
        num_last_seg_extra_const = 3; %position, velocity, accleration = 0 at end
        [A_EQ_remaining_seg_1, B_eq_remaining_seg_1] = getPerSegmentConstraints(zeros(num_seg_1_extra_const,1), 2:2+num_seg_1_extra_const-1, 1, 0);
        [A_EQ_remaining_seg_last, B_eq_remaining_seg_last] = getPerSegmentConstraints(zeros(num_last_seg_extra_const,1), 2:2+num_last_seg_extra_const-1, num_segments, time_per_segments(num_segments));
        A_EQ=[A_EQ_segment_boundary;A_EQ_position;A_EQ_remaining_seg_1;A_EQ_remaining_seg_last];
        B_eq=[B_eq_segment_boundary;B_eq_position;B_eq_remaining_seg_1;B_eq_remaining_seg_last];
        [A_INEQ, B_ineq] = get_corridorConstraints(waypoints(i,:), corridorWidths, noCorridorConstraints);
        %parameters(:,i) = A_EQ\B_eq;
        %Solving constraint optimization problem to find parameters
        %objectiveFun can be passed, but in this code it is hard coded to
        %integration(square(snap)). Integration can be done analytically in this
        %case and used. But we are numerically integrating square(snap)
        objFunction = @(params) integral(@(t) objectiveFunValue(t, snapDerivativeIndex, params), 0,time_per_segments(end));

        %Col vector for all segements [param_seg1, param_seg2...]'
        initialParams = zeros(num_segments*total_element, 1);
        options = optimoptions('fmincon','Algorithm','sqp');
        [params, fval, exitflag, output] = fmincon(objFunction, initialParams, A_INEQ, B_ineq, A_EQ, B_eq,[],[],[],options);
        disp(output);
        parameters(:,i) = params
    end
    
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
    v = zeros(1,total_length);
    v(start_pos:start_pos+length(x)-1) = x;
end

function [A_EQ, B_eq] = getContinuityConstraintForSegmentBoundary(maxSmoothness)
    global time_per_segments derivative_cache total_element num_segments;

    matlabOffset = 1;

    A_EQ = zeros((num_segments - 1)*maxSmoothness, total_element*num_segments);
    B_eq = zeros((num_segments - 1)*maxSmoothness, 1);

    for index= 1:num_segments-1
        time = time_per_segments(index);
        start = (index -1)*total_element + matlabOffset;
        for smoothDerivative =1:maxSmoothness
            constraintRow = derivative_cache{smoothDerivative}(time);

            %The way we make sure derivative is same at segment boundary is by
            %[polynomial_derivative, - polnomial_derivate,...] *
            %[params1,params2,...]' = [0...]. Thats why B_eq is all zero
            %Also we generate each row for as many smooth derivative like
            %position, velocity, acceleration etc. we want
            A_EQ((index - 1)*maxSmoothness + smoothDerivative, start:start+total_element*2-1) = [constraintRow, -constraintRow];
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

function start_index_segment = getSegStartFromConstraintChain(waypointIndex)
    global total_element
    seg_no = (waypointIndex-1);
    if waypointIndex==1
       seg_no = 1;
    end
    start_index_segment = (seg_no-1)*total_element+1;
end

%Corridor constraints
function waypointsUV = get_unitVectors_waypoints(duration_per_segment, waypoints)
    d = diff(waypoints);
    length = sqrt(d.^2 + duration_per_segment.^2);
    waypointsUV = [(d./length); (duration_per_segment./length)];
end

function [startPoint, endPoint] = get_segment(time, waypoints)
    segNo = get_segment_no(time);
    startPoint = waypoints(segNo);
    endPoint = waypoint(segNo+1);
end

function uv = get_unitVector(time, waypointsUV)
    segNo = get_segment_no(time);
    uv = waypointsUV(:, segNo);
end

function segNo = get_segment_no(time)
    global time_per_segments
    segNo = find(time_per_segments >= time,1);
end

function [coordinateWidth, posWidth, timeWidth] = get_adjusted_coordinate_width(actualWidth)
    noOfCoordinates = 3;
    coordinateWidth = actualWidth/sqrt(noOfCoordinates);
    posWidth = coordinateWidth/sqrt(2);
    timeWidth = coordinateWidth/sqrt(2);
end

function [A_INEQ, B_ineq] = get_corridorConstraints(waypoints, corridorWidths, noConstraintsPerSegment)
    global time_per_segments num_segments derivative_cache total_element
    duration_per_segments = [time_per_segments(1) diff(time_per_segments)];
    waypointsUV = get_unitVectors_waypoints(duration_per_segments, waypoints)
    A_INEQ = [];
    B_ineq = [];
    A_INEQ_field_B_ineq = [];
    position = derivative_cache{1};
    startTimeSeg = 0;
    for seg=1:num_segments
        [cW, pW, tW] = get_adjusted_coordinate_width(corridorWidths(seg));
        tStep = duration_per_segments(seg)/(noConstraintsPerSegment+1);
        for t = tStep:tStep:duration_per_segments(seg)-tStep
            constraintTime = startTimeSeg + t;
            generateCorridorConstraint();
        end
        startTimeSeg = time_per_segments(seg);
    end

    function generateCorridorConstraint()
        uv = get_unitVector(constraintTime, waypointsUV);
        u1 = uv(1);
        u2 = uv(2);
        y2 = constraintTime;
        x2 = startTimeSeg;
        x1 = waypoints(seg);
        pcc_aineq_bineq =  ((y2-x2)*u1*u2+x1*u1^2) + x1;
        tcc_aineq_bineq =  -(y2+x1*u1*u2-(y2-x2)*u2^2) + x2;
        pcc_aineq = position(constraintTime).*(1-u1^2);
        tcc_aineq = -position(constraintTime).*(u1*u2);
        seg_start_index = (seg-1)*total_element+1;
        aineq_length = num_segments*total_element;
        A_INEQ = [A_INEQ;pad_with_zero(pcc_aineq, seg_start_index, aineq_length);pad_with_zero(tcc_aineq, seg_start_index, aineq_length)];
        A_INEQ_field_B_ineq = [A_INEQ_field_B_ineq; pcc_aineq_bineq; tcc_aineq_bineq];
        B_ineq = [B_ineq;pW;tW];
    end

    %Make absolute constraints into 2 inequality constraints
    %replicate each row for > and <, since comparison is for absolute
    %values
    A_INEQ = kron(A_INEQ, [1; -1]);
    B_ineq = kron(B_ineq, [1; 1]) + kron(A_INEQ_field_B_ineq, [1;-1]); %because we shift field from A_INEQ to b_ineq, it would be 1 and -1 unlike 1,1 for b_ineq
end

function objVals = objectiveFunValue(times, objectiveDerivative, parameters)
    global derivative_cache total_element num_segments 
    A_EQ = zeros(size(times,2), total_element*num_segments);
    for index=1:size(times,2)
        t = times(index);
        segNo = get_segment_no(t);
        segStartIndex = (segNo-1)*total_element+1;
        A_EQ(index, segStartIndex:segStartIndex+total_element-1) =  derivative_cache{objectiveDerivative}(t);
    end
    objVals = (A_EQ*parameters)'.^2;
end









