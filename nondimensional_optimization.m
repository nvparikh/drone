function Parameter = nondimensional_optimization(objectiveFun, Parameter_0, tmin, tmax, pos_tmin, pos_tmax, velocity_tmin, velocity_tmax, lb, ub, maxCorridor)
   A_EQ = generate_A_EQ(tmin, tmax);
   B_eq = [pos_tmin; pos_tmax; velocity_tmax];
   
   [A_INEQ, B_ineq] = generatePositionCorridorConstraints(tmax, tmin, pos_tmax, pos_tmin, maxCorridor);
   integral_objectiveFun = @(tParameter) numericIntegration(@(t) objectiveFun(t, tParameter), tmin, tmax);
   Parameter = fmincon(integral_objectiveFun, Parameter_0, A_INEQ, B_ineq, A_EQ, B_eq, lb, ub);
   t = tmin:0.05:tmax;
   position = fliplr(Parameter) * [t.^6; t.^5; t.^4; t.^3; t.^2; t; ones(1, size(t,2))];
   figure(1);
   plot(t, position);
   title('position');
   figure(2);
   velocity = fliplr(Parameter) * [6*t.^5; 5*t.^4; 4*t.^3; 3*t.^2; 2*t; ones(1, size(t,2)); zeros(1,size(t,2))];
   plot(t, velocity);
   title('velocity');
   figure(3);
   acceleration = fliplr(Parameter) * [5*6*t.^4; 4*5*t.^3; 3*4*t.^2; 2*3*t; 2*ones(1, size(t,2)); zeros(1,size(t,2)); zeros(1,size(t,2))];
   plot(t,acceleration);
   title('acceleration');
   figure(4)
   snap_square = objectiveFun(t, Parameter);
   plot(t, snap_square);
   title('snap_square');
end

function integration = numericIntegration(fun, tmin, tmax)
    integration = integral(fun, tmin, tmax);
end

function A_EQ = generate_A_EQ(tmin, tmax)
    A_EQ = [1 tmin tmin^2 tmin^3 tmin^4 tmin^5 tmin^6; 
            1 tmax tmax^2 tmax^3 tmax^4 tmax^5 tmax^6;
            %0 1 2*tmin 3*tmin^2 4*tmin^3 5*tmin^4 6*tmin^5;
            0 1 2*tmax 3*tmax^2 4*tmax^3 5*tmax^4 6*tmax^5];
end

% A_INEQ*Parameter - 
%(A_INEQ*Parameter . UnitVector(along trajectory)) UnitVector^)
% < maxCorridor
function [A_INEQ, B_ineq] = generatePositionCorridorConstraints(tmax, tmin , pmax, pmin, maxCorridor)
    no_constraints = 30;
    unitVector = get_unitVector_straight_path(pmax, pmin, tmax, tmin);
    tConstraints = tmin:(tmax - tmin)/(no_constraints - 1):tmax;
    A_INEQ = zeros(no_constraints, 7);
    for iter = 1:no_constraints
        t = tConstraints(iter);
        A_INEQ(iter,:) = [1 t t^2 t^3 t^4 t^5 t^6];
    end
    A_INEQ = A_INEQ * (1 - unitVector(1)^2);
    %replicate each row for > and <, since comparison is for absolute
    %values
    A_INEQ = kron(A_INEQ, [1; -1]);
    
    B_ineq = (maxCorridor*ones(no_constraints, 1));
    B_ineq = kron(B_ineq, [1; 1]) + kron((tConstraints' - tmin)*unitVector(1)*unitVector(2), [1; -1]);% |pos| <= max_corridor
end

function unitVector = get_unitVector_straight_path(pmax, pmin, tmax, tmin)
    unitVector_straight_path_magnitude = sqrt((pmax - pmin)^2 + (tmax - tmin)^2);
    unitVector = [(pmax - pmin) / unitVector_straight_path_magnitude (tmax - tmin)/unitVector_straight_path_magnitude];
end

function val = test(t)
  val = [20*360 1*240 34*120 0 0 0 0] * [t.^2 t.^ t.^4 t.^3 t.^2 t 1]';
end
