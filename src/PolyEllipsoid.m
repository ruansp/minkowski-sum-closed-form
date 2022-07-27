classdef PolyEllipsoid < handle
    %PolyEllipsoid class to construct poly-ellipsoid shape, 3D.
    %
    % Author:
    %    Sipu Ruan, 2021
    
    properties
        a      % 6D array reflecting the shape
        tc     % Center point
        q      % Quaternion for rotation
        
        N      % Number of parameterized points on surface
        
        vel     % velocity of the object
    end
    
    properties
        eps = [0,0] % Indicator for ellipsoid
    end
    
    properties (Access = private)
        eta    % Meshgrid struct of eta angles
        omega  % Meshgrid struct of omega angles
    end
    
    methods
        %% Constructor
        function obj = PolyEllipsoid(val)
            if length(val) < 4
                error('Length of inputs at least 4!')
            elseif length(val{1}) ~= 6
                error('Semi-axes length size required to equal 6!')
            elseif length(val{2}) ~= 3
                error('Center point size required to equal 3!')
            elseif length(val{3}) ~= 4
                error('Quaternion size required to equal 4!')
            else
                obj.a = val{1};
                obj.tc = val{2};
                obj.q = val{3};
                
                obj.N = val{4};
                [obj.eta, obj.omega] = meshgrid(...
                    -pi/2:pi/(obj.N(1)-1):pi/2,...
                    -pi-1e-6:2*pi/(obj.N(2)-1):pi+1e-6);
                
                if length(val) == 5
                    if length(val{5}) ~= 6
                        error('Velocity must be 6x1 array: [linear; angular (in exponential coordinates)]');
                    else
                        obj.vel = val{5};
                    end
                end
            end
        end

        
        %% GetImplicitFunction implicit function
        function f = GetImplicitFunction(obj, x)
            aa = obj.GetShapeParamOctant(x);
            f = x' * diag(aa.^(-2)) * x - 1;
        end
        
        %% Get surface point cloud
        function pnt_t = GetPoints(obj)
            % GetPoints Get point cloud with pose
            
            pnt_t = quat2rotm(obj.q) * obj.GetPointsCanonical() + obj.tc;
        end
        
        function [x_t, y_t, z_t] = GetSurfPoints(obj)
            % GetSurfPoints Get surface points with pose
            
            pnt_t = obj.GetPoints();
            
            x_t = reshape(pnt_t(1,:), obj.N(1), obj.N(2));
            y_t = reshape(pnt_t(2,:), obj.N(1), obj.N(2));
            z_t = reshape(pnt_t(3,:), obj.N(1), obj.N(2));
        end
        
        function [x, y, z] = GetSurfPointsCanonical(obj)
            % GetSurfPointsCanonical Get sampled points on canonical
            %  surface
            
            % Generate unit hyperspherical coordinates
            [u_x, u_y, u_z] = obj.GetSurfHypersphere();
            
            % Decide octant
            aa = obj.a(1) * ones(size(u_x));
            aa(u_x < 0) = obj.a(2);
            bb = obj.a(3) * ones(size(u_y));
            bb(u_y < 0) = obj.a(4);
            cc = obj.a(5) * ones(size(u_z));
            cc(u_z < 0) = obj.a(6);
            
            % Canonical coordinates
            x = aa .* u_x;
            y = bb .* u_y;
            z = cc .* u_z;
        end
        
        function pnt = GetPointsCanonical(obj)
            % GetPointsCanonical Get point cloud on surface, in
            % canonical form
            
            [x, y, z] = obj.GetSurfPointsCanonical();
            
            xx = reshape(x, 1, prod(obj.N));
            yy = reshape(y, 1, prod(obj.N));
            zz = reshape(z, 1, prod(obj.N));
            
            pnt = [xx;yy;zz];
        end
        
        %% Get boundary points from different parameterizations
        function f_m = GetPointsFromGradient(obj, m)
            % GetPointsFromGradient get boundary points from un-normalized
            % gradient
            %  Input:
            %    m: gradients, 3xN vector
            
            num = size(m,2);
            f_m = nan(3,num);
            
            for i = 1:num
                aa = obj.GetShapeParamOctant(m(:,i));
                f_m(:,i) = 1/2 * diag(aa)^2 * m(:,i);
            end
        end
        
        function f_n = GetPointsFromNormal(obj, n)
            % GetPointsFromNormal Get boundary points from outward unit
            % normals
            %  Input:
            %    n: normals, 3xN vector
            
            num = size(n,2);
            f_n = nan(3,num);
            
            for i = 1:num
                aa = obj.GetShapeParamOctant(n(:,i));
                f_n(:,i) = diag(aa)^2 * n(:,i) ./ norm(diag(aa) * n(:,i));
            end
        end
        
        %% Get surface gradient cloud
        function m_t = GetGradients(obj)
            m_t = quat2rotm(obj.q) * obj.GetGradientsCanonical();
        end
        
        function [m_x, m_y, m_z] = GetSurfGradientsCanonical(obj)
            % GetSurfGradientsCanonical Get gradient vectors
            % (un-normalized) on boundary, canonical form, as surface
            
            % Generate unit hyperspherical coordinates
            [u_x, u_y, u_z] = obj.GetSurfHypersphere();
            
            % Decide octant
            aa = obj.a(1) * ones(size(u_x));
            aa(u_x < 0) = obj.a(2);
            bb = obj.a(3) * ones(size(u_y));
            bb(u_y < 0) = obj.a(4);
            cc = obj.a(5) * ones(size(u_z));
            cc(u_z < 0) = obj.a(6);
            
            m_x = 2./aa .* u_x;
            m_y = 2./bb .* u_y;
            m_z = 2./cc .* u_z;
        end
        
        function m = GetGradientsCanonical(obj)
            % GetGradientsCanonical Get gradient vectors (un-normalized) on
            %  boundary, canonical form
            
            [m_x, m_y, m_z] = obj.GetSurfGradientsCanonical();
            
            % Reshape to 3xN array
            [mm, nn] = size(m_x);
            m = nan(3,mm*nn);
            m(1,:) = reshape(m_x, 1, mm*nn);
            m(2,:) = reshape(m_y, 1, mm*nn);
            m(3,:) = reshape(m_z, 1, mm*nn);
        end
        
        %% Get gradients from different parameterizations
        function m = GetGradientsFromSpherical(obj, th)
            % GetGradientsFromCartesian Get surface gadient from cartisian
            % coordinates
            %  Input:
            %    th: array of spherical coordinates, Nx2 vector

            u(1,:) = cos(th(:,1)) .* cos(th(:,2));
            u(2,:) = cos(th(:,1)) .* sin(th(:,2));
            u(3,:) = sin(th(:,1));
            
            m = GetGradientsFromHypersphere(obj, u);
        end
        
        function m = GetGradientsFromCartesian(obj, x)
            % GetGradientsFromCartesian Get surface gadient from cartisian
            % coordinates
            %  Input:
            %    x: cartesian coordinates of points, 3xN vector
            
            num = size(x,2);
            m = nan(3,num);
            
            for i = 1:num
                aa = obj.GetShapeParamOctant(x(:,i));
                m(:,i) = 2 * diag(aa.^(-2)) * x(:,i);
            end
        end
        
        function m = GetGradientsFromHypersphere(obj, u)
            % GetGradientsFromHypersphere Get gradient on boundary based on
            % unit hypersphere parameter
            %  Input:
            %    u: unit hypersphere parameters, 3xN vector
            
            num = size(u,2);
            m = nan(3,num);
            
            for i = 1:num
                aa = obj.GetShapeParamOctant(u(:,i));
                m(:,i) = 2 * diag(aa.^(-1)) * u(:,i);
            end
        end
        
        function m = GetGradientsFromDirection(obj, d)
            % GetGradientsFromHypersphere Get gradient on boundary based on
            % unit hypersphere parameter
            %  Input:
            %    d: directional vector (normalization not required),
            %       3xN vector
            
            n = d./sqrt(sum(d.^2,1));
            num = size(n,2);
            m = nan(3,num);
            
            for i = 1:num
                aa = obj.GetShapeParamOctant(n(:,i));
                m(:,i) = 2 * n(:,i) / norm(diag(aa) * n(:,i));
            end
        end   
        
        %% Get hypersphere parameters
        function [u_x, u_y, u_z] = GetSurfHypersphere(obj)
            u_x = cos(obj.eta) .* cos(obj.omega);
            u_y = cos(obj.eta) .* sin(obj.omega);
            u_z = sin(obj.eta) .* ones(size(obj.omega));
        end
        
        function u = GetHypersphereFromGradient(obj, m)
            % GetHypersphereFromGradient Get unit hypersphere parameters
            % from gradient
            %  Input:
            %    m: gradients, 3xN vector
            
            num = size(m,2);
            u = nan(3,num);
            
            for i = 1:num
                aa = obj.GetShapeParamOctant(m(:,i));
                u(:,i) = 1/2 * diag(aa) * m(:,i);
            end
        end
        
        %% GetSurf Get mesh information as surface object
        function S = GetSurf(obj)
            [x, y, z] = obj.GetSurfPoints();
            
            S = surf(x, y, z, 'FaceColor', 'none', 'EdgeColor', 'none');
        end
        
        % Plot surface
        function PlotShape(obj, color, faceAlpha)
            [x, y, z] = obj.GetSurfPoints();
            
            surf(x, y, z, 'FaceColor', color, 'EdgeColor', 'none',...
                'FaceAlpha', faceAlpha);
        end
        
        %% Get ellipsoid for the corresponding octant info
        function aa = GetShapeParamOctant(obj, pnt)
            aa = [obj.a(1), obj.a(3), obj.a(5)];
            if pnt(1) < 0; aa(1) = obj.a(2); end
            if pnt(2) < 0; aa(2) = obj.a(4); end
            if pnt(3) < 0; aa(3) = obj.a(6); end
        end
    end
end