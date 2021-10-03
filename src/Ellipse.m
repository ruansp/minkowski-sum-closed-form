classdef Ellipse < SuperEllipse & EllipsoidND
    % Ellipse class, 2D Ellipse. Derived from EllipsoidND and SuperEllipse 
    % classes
    %
    %  Author:
    %    Sipu Ruan, ruansp@nus.edu.sg, 2021
    %
    %  See also
    %    SuperEllipse, EllipsoidND

    methods
        %% Constructor
        function obj = Ellipse(val)
            obj@SuperEllipse([val(1:2), 1, val(3:end)]);
            obj@EllipsoidND(diag(val(1:2)), rot2(val(6)));
        end
        
        %% GetImplicitFunction implicit function
        function f = GetImplicitFunction(obj, x)
            f = obj.GetImplicitFunction@EllipsoidND(x);
        end
        
        %% Get boundary points from different parameterizations
        function f_m = GetPointsFromGradient(obj, m)
            % GetPointsFromGradient get boundary points from un-normalized
            % gradient
            %  Input:
            %    m: gradients, 2xN vector
            
            f_m = obj.GetPointsFromGradient@EllipsoidND(m);
        end
        
        function f_n = GetPointsFromNormal(obj, n)
            % GetPointsFromNormal Get boundary points from outward unit
            % normals
            %  Input:
            %    n: normals, 2xN vector
            
            f_n = obj.GetPointsFromNormal@EllipsoidND(n);
        end
        
        %% Get gradients from different parameterizations
        function m = GetGradientsFromCartesian(obj, x)
            % GetGradientsFromCartesian Get surface gadient from cartisian
            % coordinates
            %  Input:
            %    x: cartesian coordinates of points, 2xN vector
            
            m = obj.GetGradientsFromCartesian@EllipsoidND(x);
        end
        
        function m = GetGradientsFromHypersphere(obj, u)
            % GetGradientsFromHypersphere Get gradient on boundary based on
            % unit hypersphere parameter
            %  Input:
            %    u: unit hypersphere parameters, 2xN vector
            
            m = obj.GetGradientsFromHypersphere@EllipsoidND(u);
        end       
        
        %% Get hypersphere parameters
        function u = GetHypersphereFromGradient(obj, m)
            % GetHypersphereFromGradient Get unit hypersphere parameters
            % from gradient
            %  Input:
            %    m: gradients, 2xN vector
            
            u = obj.GetHypersphereFromGradient@EllipsoidND(m);
        end
    end
end