classdef EllipsoidND < handle
    % EllipsoidND class, n-dimensional ellipsoid with different operations.
    %
    %  Inputs:
    %    VA: Shape matrix in canonical form, diagonal matrix with semi-axis
    %        lengths in diagonal entries
    %    R : Rotation matrix of the ellipsoid
    %
    %  Author:
    %    Sipu Ruan, ruansp@nus.edu.sg, 2021
    
    properties
        dim % Dimension
        VA  % Shape matrix in canonical form
        A   % Shape matrix in rotated form
    end
    
    methods
        %% Constructor
        function obj = EllipsoidND(VA, R)
            obj.VA = VA;
            obj.A = R * obj.VA * R';
        end
        
        %% GetImplicitFunction implicit function
        function f = GetImplicitFunction(obj, x)
            % GetImplicitFunction Get the implicit expression value of a
            % point in global frame
            %
            %  Inputs:
            %    x: A set of points in space, nxN vector
            %
            %  Outputs:
            %    f: Implicit expressions of ellipsoid, 1xN vector
            
            f = x' * diag( obj.a.^(-2) ) * x - 1;
        end
        
        %% Get boundary points from different parameterizations
        function f_m = GetPointsFromGradient(obj, m)
            % GetPointsFromGradient get boundary points from un-normalized
            % gradient
            %
            %  Inputs:
            %    m  : Gradients, nxN vector
            %
            %  Outputs:
            %    f_m: Surface points according to the gradients, nxN vector
            
            f_m = 1/2 * obj.VA.^2 * m;
        end
        
        function f_n = GetPointsFromNormal(obj, n)
            % GetPointsFromNormal Get boundary points from outward unit
            % normals
            %
            %  Inputs:
            %    n  : Normals, nxN vector
            % 
            %  Outputs:
            %    f_n: Surface points corresponding to the normals, nxN
            %         vector
            
            f_n = obj.VA.^2 * n ./ sqrt( sum((obj.VA * n).^2, 1) );
        end
        
        %% Get gradients from different parameterizations
        function m = GetGradientsFromCartesian(obj, x)
            % GetGradientsFromCartesian Get surface gadient from cartisian
            % coordinates
            %
            %  Inputs:
            %    x: Cartesian coordinates of points, nxN vector
            %
            %  Outputs:
            %    m: Gradient vectors of the points, nxN vector
            
            m = diag( 2./obj.a.^2 ) * x;
        end
        
        function m = GetGradientsFromHypersphere(obj, u)
            % GetGradientsFromHypersphere Get gradient on boundary based on
            % unit hypersphere parameter
            %
            %  Inputs:
            %    u: Unit hypersphere parameters, nxN vector
            %
            %  Outputs:
            %    m: Gradients, nxN vector
            
            m = diag( 2./obj.a ) * u;
        end       
        
        %% Get hypersphere parameters
        function u = GetHypersphereFromGradient(obj, m)
            % GetHypersphereFromGradient Get unit hypersphere parameters
            % from gradient
            %
            %  Inputs:
            %    m: Gradients, nxN vector
            %
            %  Outputs:
            %    u: Unit hypersphere parameters, nxN vector
            
            u = 1/2 * obj.VA * m;
        end
    end
end