classdef MinkSumClosedForm < handle
    % MinkSumClosedForm Computes exact closed-form Minkowski sums between 
    % two convex bodies with smooth and positively curved boundaries
    %
    %  Inputs:
    %    s1, s2: Geometry class
    %    M1, M2: Linear transformation matrix
    %
    %  Author:
    %    Sipu Ruan, ruansp@nus.edu.sg, 2021
    %
    %  See also
    %    Ellipse, Ellipsoid, SuperEllipse, SuperQuadrics
    
    %% Variables
    properties
        s1    % Class of convex smooth body
        s2    % Class of convex smooth body
        M1    % Linear transformations for s1
        M2    % Linear transformations for s2
        dim   % Dimension of space
    end
    
    properties (Access = private)
        is_s2_e = false  % Whether s2 is ellipsoid
    end
    
    %% Functions
    methods
        %% Constructor
        function obj = MinkSumClosedForm(s1, s2, M1, M2)
            obj.s1 = s1;
            obj.s2 = s2;
            obj.M1 = M1;
            obj.M2 = M2;
            obj.dim = size(obj.M1,1);
            
            if ( obj.dim == 2 )
                if( obj.s2.eps == 1 )
                    obj.is_s2_e = true;
                    obj.s2 = Ellipse([s2.a, 0, s2.tc', s2.ang, s2.N]);
                end
            elseif( obj.dim == 3 )
                if( obj.s2.eps(1) == 1 && obj.s2.eps(2) == 1 )
                    obj.is_s2_e = true;
                    obj.s2 = Ellipsoid({s2.a, [0,0], s2.tc', s2.q, s2.N});
                end
            end
        end
        
        %% GetMinkSumFromGradient compute Minkowski sums from gradients
        % Inputs:
        %   m1  : gradient vectors on s1
        %
        % Output:
        %   mink: boundary points of closed-form Minkowski sum
        function mink = GetMinkSumFromGradient(obj, m1)
            % Compute f_1(m_1)
            f_m1 = obj.s1.GetPointsFromGradient( m1 );
            
            % Linear deformations
            m1_t = obj.M2' * (obj.M1' \ m1);
            
            % Compute f_2(m_1)
            if obj.is_s2_e
                % s2 is ellipsoid
                A2 = obj.s2.VA;
                f_m2 = A2^2 * (-m1_t) ./ sqrt( sum((A2*m1_t).^2, 1) );
            else
                % s2 is not ellipsoid
                % Phi(m1)
                phi_m1_t = obj.Phi(-m1_t);
                
                f_m2 = obj.s2.GetPointsFromGradient( ...
                    - phi_m1_t ./ sqrt( sum(m1_t.^2, 1) ) .* m1_t );
            end
            
            % Compute Minkowski sum
            mink = obj.M1 * f_m1 - obj.M2 * f_m2;
        end
        
        % Phi |m2| = Phi(m1) Compute |m2| as a function of m1
        function phi_m1 = Phi(obj, m1)
            g2_m1_inv = obj.s2.GetHypersphereFromGradient(m1);
            m2 = obj.s2.GetGradientsFromHypersphere( ...
                g2_m1_inv ./ sqrt( sum(g2_m1_inv.^2, 1) ) );
            
            phi_m1 = sqrt( sum(m2.^2, 1) );
        end
        
        %% GetMinkSumFromNormal compute Minkowski sums from normals
        % Inputs:
        %   n1  : normal vectors on s1
        %
        % Output:
        %   mink: boundary points of closed-form Minkowski sum
        function mink = GetMinkSumFromNormal(obj, n1)
            % Compute f_i(n_1)
            f_n1 = obj.s1.GetPointsFromNormal( n1 );
            f_n2 = obj.s2.GetPointsFromNormal( obj.M2' * (obj.M1'\(-n1)) );
            
            mink = obj.M1 * f_n1 - obj.M2 * f_n2;
        end
        
        %% GetMinkSumGeometric compute Minkowski sums from geometric method
        % Output:
        %   mink: boundary points of closed-form Minkowski sum
        %
        % Note:
        %   Works only when one object is ellipsoid
        function mink = GetMinkSumGeometric(obj)
            if ~obj.is_s2_e
                error('Second body should be an ellipsoid!');
            end
            
            % Shrinking transformation
            Tinv = obj.M2 * obj.s2.VA / obj.M2;
            
            % Gradient of s1
            gradPhi = obj.s1.GetGradients();
            gradPhi_t = obj.M2 * obj.M2' * (obj.M1' \ gradPhi);
            
            % Closed-Form Minkowski Sum
            mink = obj.M1 * obj.s1.GetPoints() + Tinv^2 * gradPhi_t ./...
                sqrt(sum( (obj.M2 \ (Tinv * gradPhi_t)).^2, 1 ));
        end
    end
end
