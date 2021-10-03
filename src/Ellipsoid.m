classdef Ellipsoid < SuperQuadrics & EllipsoidND
    % Ellipsoid class, 3D ellipsoid. Derived from EllipsoidND and
    % SuperQuadrics classes
    %
    %  Inputs:
    %    val{1}: Semi-axis length, 1x3 vector
    %    val{2}: Tapering factor, 1x2 vector
    %    val{3}: Center coordinate, 3x1 vector
    %    val{4}: Quaternion for orientation, 1x4 vector (normalized)
    %    val{5}: Indicators of number of points on surface, 1x2 vector
    %
    %  Author:
    %    Sipu Ruan, ruansp@nus.edu.sg, 2021
    %
    %  See also
    %    SuperQuadrics, EllipsoidND

    methods
        %% Constructor
        function obj = Ellipsoid(val)
            obj@SuperQuadrics([val{1}, [1,1], val(2:end)]);
            obj@EllipsoidND(diag(val{1}), quat2rotm(val{4}));
        end
        
        %% GetImplicitFunction implicit function
        function f = GetImplicitFunction(obj, x)
            % GetImplicitFunction Get implicit expression of a point in
            % space
            %
            %  Inputs:
            %    x: A set of points in space, 3xN vector
            %
            %  Outputs:
            %    f: Implicit expressions of the point set, 1xN vector
            
            f = obj.GetImplicitFunction@EllipsoidND(x);
        end
        
        %% Get boundary points from different parameterizations
        function f_m = GetPointsFromGradient(obj, m)
            % GetPointsFromGradient get boundary points from un-normalized
            % gradient
            %
            %  Inputs:
            %    m  : gradients, 3xN vector
            %
            %  Outputs:
            %    f_m: Surface points, 3xN vector
            
            f_m = obj.GetPointsFromGradient@EllipsoidND(m);
        end
        
        function f_n = GetPointsFromNormal(obj, n)
            % GetPointsFromNormal Get boundary points from outward unit
            % normals
            % 
            %  Inputs:
            %    n: Normals, 3xN vector
            % 
            %  Outputs:
            %    f_n: Surface points, 3xN vector
            
            f_n = obj.GetPointsFromNormal@EllipsoidND(n);
        end
        
        %% Get gradients from different parameterizations
        function m = GetGradientsFromCartesian(obj, x)
            % GetGradientsFromCartesian Get surface gadient from cartisian
            % coordinates
            %
            %  Inputs:
            %    x: Cartesian coordinates of points, 3xN vector
            %
            %  Outputs:
            %    m: Gradients of the points, 3xN vector
            
            m = obj.GetGradientsFromCartesian@EllipsoidND(x);
        end
        
        function m = GetGradientsFromHypersphere(obj, u)
            % GetGradientsFromHypersphere Get gradient on boundary based on
            % unit hypersphere parameter
            %
            %  Inputs:
            %    u: Unit hypersphere parameters, 3xN vector
            %
            %  Outputs:
            %    m: Gradients, 3xN vector
            
            m = obj.GetGradientsFromHypersphere@EllipsoidND(u);
        end       
        
        %% Get hypersphere parameters
        function u = GetHypersphereFromGradient(obj, m)
            % GetHypersphereFromGradient Get unit hypersphere parameters
            % from gradient
            %
            %  Inputs:
            %    m: Gradients, 3xN vector
            %
            %  Outputs:
            %    u: Unit hypersphere parameters, 3xN vector
            
            u = obj.GetHypersphereFromGradient@EllipsoidND(m);
        end
    end
end