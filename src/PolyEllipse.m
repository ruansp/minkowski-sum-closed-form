classdef PolyEllipse < handle
    %PolyEllipse class to construct poly-ellipse shape, 2D.
    %
    % Author:
    %    Sipu Ruan, 2021
    
    properties
        a       % 4D array reflecting the shape
        tc      % Center point
        ang     % Angle of rotation
        
        N       % Number of parameterized points on surface
    end
    
    properties
        eps = 0 % Indicator for ellipse
    end
    
    properties (Access = private)
        th      % theta angles
    end
    
    methods
        %% Constructor
        function obj = PolyEllipse(val)
            if length(val) < 4
                error('Length of inputs at least 4!')
            else
                if length(val{1}) ~= 4
                    error('Semi-axes length size required to equal 4!')
                elseif length(val{2}) ~= 2
                    error('Center point size required to equal 2!')
                elseif length(val{3}) ~= 1
                    error('Angle size required to equal 1!')
                else
                    obj.a = val{1};
                    obj.tc = val{2};
                    obj.ang = val{3};

                    obj.N = val{4};
                    obj.th = -pi-1e-6:2*pi/(obj.N-1):pi+1e-6;
                end
            end
        end
        
        %% GetImplicitFunction implicit function
        function f = GetImplicitFunction(obj, x)
            aa = obj.GetShapeParamQuadrant(x);
            f = x' * diag(aa.^(-2)) * x - 1;
        end
        
        %% Get surface point cloud
        function x_t = GetPoints(obj)
            % GetPoints Get point cloud with pose
            
            x_t = rot2(obj.ang) * obj.GetPointsCanonical() + obj.tc;
        end
        
        function x = GetPointsCanonical(obj)
            % GetPointsCanonical Get point cloud on surface, in
            % canonical form
            
            % Generate unit hyperspherical coordinates
            u = obj.GetHypersphere();
            
            % Decide quadrant
            aa = obj.a(1) * ones(1, size(u,2));
            aa(u(1,:) < 0) = obj.a(2);
            bb = obj.a(3) * ones(1, size(u,2));
            bb(u(2,:) < 0) = obj.a(4);
            
            % Gradient
            x = nan(2,obj.N);
            x(1,:) = aa .* u(1,:);
            x(2,:) = bb .* u(2,:);
        end
        
        %% Get boundary points from different parameterizations
        function f_m = GetPointsFromGradient(obj, m)
            % GetPointsFromGradient get boundary points from un-normalized
            % gradient
            %  Input:
            %    m: gradients, 2xN vector
            
            num = size(m,2);
            f_m = nan(2,num);
            
            for i = 1:num
                aa = obj.GetShapeParamQuadrant(m(:,i));
                f_m(:,i) = 1/2 * diag(aa)^2 * m(:,i);
            end
        end
        
        function f_n = GetPointsFromNormal(obj, n)
            % GetPointsFromNormal Get boundary points from outward unit
            % normals
            %  Input:
            %    n: normals, 2xN vector
            
            num = size(n,2);
            f_n = nan(2,num);
            
            for i = 1:num
                aa = obj.GetShapeParamQuadrant(n(:,i));
                f_n(:,i) = diag(aa)^2 * n(:,i) ./ norm(diag(aa) * n(:,i));
            end
        end
        
        %% Get surface gradient cloud
        function m_t = GetGradients(obj)
            m_t = rot2(obj.ang) * obj.GetGradientsCanonical();
        end
        
        function m = GetGradientsCanonical(obj)
            % GetGradientsCanonical Get gradient vectors (un-normalized) on
            %  boundary, canonical form
            
            % Generate unit hyperspherical coordinates
            u = obj.GetHypersphere();
            
            % Decide quadrant
            aa = obj.a(1) * ones(1, size(u,2));
            aa(u(1,:) < 0) = obj.a(2);
            bb = obj.a(3) * ones(1, size(u,2));
            bb(u(2,:) < 0) = obj.a(4);
            
            % Gradient
            m = nan(2,obj.N);
            m(1,:) = 2./aa .* u(1,:);
            m(2,:) = 2./bb .* u(2,:);
        end
        
        %% Get gradients from different parameterizations
        function m = GetGradientsFromCartesian(obj, x)
            % GetGradientsFromCartesian Get surface gadient from cartisian
            % coordinates
            %  Input:
            %    x: cartesian coordinates of points, 2xN vector
            
            num = size(x,2);
            m = nan(2,num);
            
            for i = 1:num
                aa = obj.GetShapeParamQuadrant(x(:,i));
                m(:,i) = 2 * diag(aa.^(-2)) * x(:,i);
            end
        end
        
        function m = GetGradientsFromHypersphere(obj, u)
            % GetGradientsFromHypersphere Get gradient on boundary based on
            % unit hypersphere parameter
            %  Input:
            %    u: unit hypersphere parameters, 2xN vector
            
            num = size(u,2);
            m = nan(2,num);
            
            for i = 1:num
                aa = obj.GetShapeParamQuadrant(u(:,i));
                m(:,i) = 2 * diag(aa.^(-1)) * u(:,i);
            end
        end       
        
        %% Get hypersphere parameters
        function u = GetHypersphere(obj)
            u(1,:) = cos(obj.th);
            u(2,:) = sin(obj.th);
        end
        
        function u = GetHypersphereFromGradient(obj, m)
            % GetHypersphereFromGradient Get unit hypersphere parameters
            % from gradient
            %  Input:
            %    m: gradients, 2xN vector
            
            num = size(m,2);
            u = nan(2,num);
            
            for i = 1:num
                aa = obj.GetShapeParamQuadrant(m(:,i));
                u(:,i) = 1/2 * diag(aa) * m(:,i);
            end
        end
        
        %% Get ellipse for the corresponding quadrant info
        function aa = GetShapeParamQuadrant(obj, pnt)
            aa = [obj.a(1), obj.a(3)];
            if pnt(1) < 0; aa(1) = obj.a(2); end
            if pnt(2) < 0; aa(2) = obj.a(4); end
        end
    end
end