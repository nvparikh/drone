%About this code
%I have accomplished generating multi path min snap trajectory
%Using MATLAB's symbolic calculations allowing only positional constraints
%Advantage of using symbolic toolbox is that it greatly reduces size of
%code however it is very very slow. So converted symbolic function into
%MATLAB annonymous function using matlabFunction which greatly increases
%spped

%Input:
%numSegment = no. of piece wise segments
%timeForEachSegment = absolute time allowed to finish each segment
%maxPositionPolynomial = what polynomial size to fit
%waypointConstraints = any positional constraints along the way
%             format = col array of [time, position],
%              0 >= time <= numSegment*timeForEachSegment
%Sample input
%waypointConstraints = [[0 0];[100 1];[0 2]];
%nondimensional_optimization(3, 1, 6, waypointConstraints)

function Parameter = nondimensional_optimization(numSegment, timeForEachSegment, maxPositionPolynomial, waypointConstraints)
global polynomialDegree
polynomialDegree = maxPositionPolynomial;

%Storing derivatives of polynomial for faster processing
global derivativeCache
derivativeCache = cell(5,1);
for index=1:5
    derivativeCache{index} = getPositionDerivatives(index-1);
end

%Make sure position, velocity, acceleration and snap remain same at segment
%boundary
[A_EQ_segmentContinum, B_eq_segmentContinum] = getContinuityConstraintForSegmentBoundary(numSegment, timeForEachSegment, 4);

%Generate constraint matrix for all positional constraints
%And combine it with smoothness constraint at segment boundary
A_EQ = [getConstraintsMatrix(waypointConstraints(:,2), 0, timeForEachSegment, numSegment); A_EQ_segmentContinum];
B_eq = [waypointConstraints(:,1); B_eq_segmentContinum];

%objectiveFun can be passed, but in this code it is hard coded to
%integration(square(snap)). Integration can be done analytically in this
%case and used. But we are numerically integrating square(snap)
integral_objectiveFun = @(tParameter) numericIntegration(@(t) objectiveFun(t, tParameter, numSegment, timeForEachSegment), 0, timeForEachSegment*numSegment);

%Col vector for all segements [param_seg1, param_seg2...]'
initialParameterValue = zeros(numSegment*(polynomialDegree + 1), 1);
Parameter = fmincon(integral_objectiveFun, initialParameterValue, [], [], A_EQ, B_eq);

%Display trajectory found
disp(Parameter);

t = 0:0.05:timeForEachSegment*numSegment;
position = getConstraintsMatrix(t', 0, timeForEachSegment, numSegment) * Parameter;
figure(1);
plot(t, position');
end

function objective = objectiveFun(t, parameter, numSegment, timeForEachSegment)
A_EQ = getConstraintsMatrix(t', 4, timeForEachSegment, numSegment);
objective = (A_EQ*parameter)'.^2;
end

function [A_EQ, B_eq] = getContinuityConstraintForSegmentBoundary(numSegment, timeForEachSegment, maxSmoothDerivativeAtBoundary)
global polynomialDegree
totalElement = polynomialDegree + 1;
global derivativeCache
matlab_offset = 1;

A_EQ = zeros((numSegment - 1)*maxSmoothDerivativeAtBoundary, totalElement*numSegment);
B_eq = zeros((numSegment - 1)*maxSmoothDerivativeAtBoundary, 1);

times = timeForEachSegment:timeForEachSegment:timeForEachSegment*numSegment;
for index= 1:size(times, 2) - 1
    time = times(index);
    start = (index -1)*totalElement + matlab_offset;
    for smoothDerivative =1:maxSmoothDerivativeAtBoundary
        constraintRow = derivativeCache{smoothDerivative}(time);
        
        %The way we make sure derivative is same at segment boundary is by
        %[polynomial_derivative, - polnomial_derivate,...] *
        %[params1,params2,...]' = [0...]. Thats why B_eq is all zero
        %Also we generate each row for as many smooth derivative like
        %position, velocity, acceleration etc. we want
        A_EQ((index - 1)*maxSmoothDerivativeAtBoundary + smoothDerivative, start:start+totalElement*2-1) = [constraintRow, -constraintRow];
    end
end
end

%Input:
%times = generate constraints for each time in times
function A_EQ_constraints = getConstraintsMatrix(times, whichDerivative, timeForEachSegment, numSegment)
global polynomialDegree
totalElement = polynomialDegree + 1;
global derivativeCache
A_EQ_constraints = zeros(size(times,1), totalElement*numSegment);
for index = 1:size(times, 1)
    time = times(index);
    segmentNo = max(0, ceil(time/timeForEachSegment) - 1); %segmentNo starts from 0, handling special case of 0
    A_EQ_constraints(index, :) =  generateRowForSegment(derivativeCache{whichDerivative+1}(time), segmentNo, numSegment);
end
end

%actualRow = row containing constraints for segments
%segmentNo = segment for which to put actualRow in.
function row = generateRowForSegment(actualRow, segmentNo, numSegment)
global polynomialDegree
totalElement = polynomialDegree + 1; %including constant
row = zeros(1, totalElement * numSegment);
matlab_offset = 1;
start = matlab_offset+segmentNo*totalElement;
row(start:start+polynomialDegree) = row(start:start + polynomialDegree) + actualRow;
end

function integration = numericIntegration(fun, tmin, tmax)
integration = integral(fun, tmin, tmax);
end

function dpositionWoParamsdx = getPositionDerivatives(numDerivatives)
global polynomialDegree
totalElement = polynomialDegree + 1; %including constant
x = [];
syms x;
positionWOParams = sym('x', [1  totalElement]);
for degree=0:polynomialDegree
    positionWOParams(degree + 1) = power(x, degree);
end
dpositionWoParamsdx = matlabFunction(diff(positionWOParams, numDerivatives));
%matlabFunction greatly improves speed
end
