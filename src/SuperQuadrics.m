classdef SuperQuadrics < handle
    % SuperQuadrics class, with different methods of operations
    %
    %  Inputs:
    %    val{1}: Semi-axis length, 1x3 vector
    %    val{2}: Exponents, 1x2 vector
    %    val{3}: Tapering factor, 1x2 vector
    %    val{4}: Center coordinate, 3x1 vector
    %    val{5}: Quaternion for orientation, 1x4 vector (normalized)
    %    val{6}: Indicators of number of points on surface, 1x2 vector
    %    val{7}: [Optional] Velocity of the object, 6x1 vector including
    %            linear and angular (in exponential coordinates) velocity
    %
    %  Author:
    %    Sipu Ruan, ruansp@nus.edu.sg, 2021
    
    properties
        a       % semi-axes lengths, 3x1 array
        eps     % exponent for the signed power function, 2x1 array
        taper   % tapering factor, 2x1 array
        
        q       % exponential coordinates for orientation, 3x1 array
        tc      % center, 3x1 array
        
        vel     % velocity of the object
        
        N       % number of points for two spherical parameters, 1x2 array
        omega   % parameter for point cloud
        eta     % parameter for point cloud
    end
    
    methods
        %% Constructor
        function obj = SuperQuadrics(val)
            if length(val) < 6
                error('Input required to be larger or equal to 6!')
            elseif length(val{1}) ~= 3
                error('Semi-axes length size required to equal 3!')
            elseif length(val{2}) ~= 2
                error('Exponents size required to equal 2!')
            elseif length(val{3}) ~= 2
                error('Tapering factor size required to equal 2!')
            elseif length(val{4}) ~= 3
                error('Center point size required to equal 3!')
            elseif length(val{5}) ~= 4
                error('Quaternion size required to equal 4!')
            else
                obj.a     = val{1};
                obj.eps   = val{2};
                obj.taper = val{3};
                
                if size(val{4},1) == 1
                    obj.tc = val{4}';
                else
                    obj.tc = val{4};
                end
                
                obj.q     = val{5};
                
                obj.N     = val{6};
                
                [obj.eta, obj.omega] = meshgrid(...
                    -pi/2:pi/(obj.N(1)-1):pi/2,...
                    -pi-1e-6:2*pi/(obj.N(2)-1):pi+1e-6);
                
                if length(val) == 7
                    if length(val{7}) ~= 6
                        error('Velocity must be 6x1 array: [linear; angular (in exponential coordinates)]');
                    else
                        obj.vel = val{7};
                    end
                end
            end
        end
        
        %% Get implicit function
        function f = GetImplicitFunction(obj, x)
            % GetImplicitFunction Get implicit function of superquadrics
            %
            %  Inputs:
            %    x: Cartesian coordinates, 3xN vector
            %
            %  Outputs:
            %    f: Implicit function \Phi(x) = 0, 1xN vector
            
            % De-Tapering deformation
            x = obj.DeTaperPoint(x);
            
            % Implicit function
            f = obj.GetImplicitFunctionCanonical(x);
        end
        
        function f = GetImplicitFunctionCanonical(obj, x)
            % GetImplicitFunctionCanonical Get implicit function of 
            %  superquadrics, in canonical form (body frame)
            %
            %  Inputs:
            %    x: Cartesian coordinates, 3xN vector
            %
            %  Outputs:
            %    f: Implicit function \Phi(x) = 0, 1xN vector
            
            % Implicit function
            f = eps_fun( ...
                eps_fun( (x(1,:) ./ obj.a(1)).^2, 1/obj.eps(2) ) +...
                eps_fun( (x(2,:) ./ obj.a(2)).^2, 1/obj.eps(2) ),...
                obj.eps(2)/obj.eps(1) ) +...
                eps_fun( (x(3,:) ./ obj.a(3)).^2, 1/obj.eps(1) ) - 1;
        end
        
        %% Get surface point cloud
        function pnt_t = GetPoints(obj)
            % GetPoints Get point cloud in global frame
            %
            %  Outputs:
            %    pnt_t: Point cloud on surface, 3xN vector
            
            pnt_can = obj.GetPointsCanonical();
            pnt = obj.TaperPoint(pnt_can);
            pnt_t = quat2rotm(obj.q) * pnt + obj.tc;
        end
        
        function [x_t, y_t, z_t] = GetSurfPoints(obj)
            % GetSurfPoints Get surface points in global frame, in the
            % structure of mesh grid
            %
            %  Outputs:
            %    x_t, y_t, z_t: Surface points, in the form of mesh grid
            
            pnt_t = obj.GetPoints();
            
            x_t = reshape(pnt_t(1,:), obj.N(1), obj.N(2));
            y_t = reshape(pnt_t(2,:), obj.N(1), obj.N(2));
            z_t = reshape(pnt_t(3,:), obj.N(1), obj.N(2));
        end
        
        function [x, y, z] = GetSurfPointsTapered(obj)
            % GetSurfPointsTapered Get sampled points on superquadric 
            % surface with tapering deformation
            %
            %  Outputs:
            %    x, y, z: Tapered surface points, in the form of mesh grid
            
            [x, y, z] = obj.GetSurfPointsCanonical();
            
            % Tapering deformation
            [x, y, z] = obj.TaperSurfPoints(x, y, z);
        end
        
        function [x, y, z] = GetSurfPointsCanonical(obj)
            % GetSurfPointsCanonical Get sampled points on canonical
            %  superquadric surface
            %
            %  Outputs:
            %    x, y, z: Surface points, in the form of mesh grid
            
            % Canonical SQ coordinates
            x = obj.a(1) .* eps_fun(cos(obj.eta), obj.eps(1))...
                .* eps_fun(cos(obj.omega), obj.eps(2));
            y = obj.a(2) .* eps_fun(cos(obj.eta), obj.eps(1))...
                .* eps_fun(sin(obj.omega), obj.eps(2));
            z = obj.a(3) .* eps_fun(sin(obj.eta), obj.eps(1))...
                .* ones(size(obj.omega));
        end
        
        function pnt = GetPointsCanonical(obj)
            % GetPointsCanonical Get point cloud of the surface, in
            % canonical form
            %
            %  Outputs:
            %    pnt: Point cloud on surface, 3xN vector
            
            [x, y, z] = obj.GetSurfPointsCanonical();
            
            xx = reshape(x, 1, prod(obj.N));
            yy = reshape(y, 1, prod(obj.N));
            zz = reshape(z, 1, prod(obj.N));
            
            pnt = [xx;yy;zz];
        end
        
        %% Get boundary points from different parameterizations
        function x = GetPointsFromParamTapered(obj, par, opt)
            % GetPointsFromParamTapered Get surface points from different
            % parameterizations
            %
            %  Inputs:
            %    par: Parameters
            %    opt: Indicator of parameterization
            %
            %  Outputs:
            %    x  : Surface points, 3xN vector
            
            if strcmp(opt, 'spherical')
                x = obj.GetPointsFromSpherical(par);
                
            elseif strcmp(opt, 'gradient')
                x = obj.GetPointsFromGradient(par);
                
            elseif strcmp(opt, 'normal')
                x = obj.GetPointsFromNormal(par);
                
            elseif strcmp(opt, 'tan-half-angle')
                x = obj.GetPointsFromTanHalfAngle(par);

            else
                error('No such parameterization for surface points!')
            end
            
            % Tapering deformation
            x = obj.TaperPoint(x);
        end
        
        function x = GetPointsFromSpherical(obj, theta)
            % GetPointsFromSpherical get boundary points from spherical
            % coordinates, canonical form
            %
            %  Inputs:
            %    theta: Angular parameters, Nx2 vector
            %
            %  Outputs:
            %    x    : Surface points, 3xN vector

            x = nan(3,size(theta,1));
            
            x(1,:) = obj.a(1) .* eps_fun(cos(theta(:,1)'), obj.eps(1))...
                .* eps_fun(cos(theta(:,2)'), obj.eps(2));
            x(2,:) = obj.a(2) .* eps_fun(cos(theta(:,1)'), obj.eps(1))...
                .* eps_fun(sin(theta(:,2)'), obj.eps(2));
            x(3,:) = obj.a(3) .* eps_fun(sin(theta(:,1)'), obj.eps(1));
        end
        
        function f_m = GetPointsFromGradient(obj, m)
            % GetPointsFromGradient get boundary points from un-normalized
            % gradient, canonical form
            %
            %  Inputs:
            %    m  : Gradients, 3xN vector
            % 
            %  Outputs:
            %    f_m: Surface point, 3xN vector
            
            f_m = nan(size(m));
            
            f_m(1,:) = obj.a(1) * eps_fun( obj.a(1)*obj.eps(1)/2 .*...
                m(1,:), obj.eps(2)/(2-obj.eps(2)) ) .*...
                eps_fun( 1 - eps_fun( (obj.a(3)*obj.eps(1)/2 * m(3,:)).^2,...
                1/(2-obj.eps(1)) ), ...
                (obj.eps(1)-obj.eps(2))/(2-obj.eps(2)) );
            f_m(2,:) = obj.a(2) * eps_fun( obj.a(2)*obj.eps(1)/2 .*...
                m(2,:), obj.eps(2)/(2-obj.eps(2)) ) .*...
                eps_fun( 1 - eps_fun( (obj.a(3)*obj.eps(1)/2 * m(3,:)).^2,...
                1/(2-obj.eps(1)) ), ...
                (obj.eps(1)-obj.eps(2))/(2-obj.eps(2)) );
            f_m(3,:) = obj.a(3) * eps_fun( obj.a(3)*obj.eps(1)/2 .*...
                m(3,:), obj.eps(1)/(2-obj.eps(1)));
        end
        
        function f_n = GetPointsFromNormal(obj, n)
            % GetPointsFromNormal Get boundary points from outward unit
            % normals, canonical form
            %
            %  Inputs:
            %    n  : Normals, 3xN vector
            %
            %  Outputs:
            %    f_n: Surface point, 3xN vector
            
            m = obj.GetGradientsFromDirection(n);
            f_n = obj.GetPointsFromGradient(m);
        end
        
        function x = GetPointsFromTanHalfAngle(obj, t)
            % GetPointsFromTanHalfAngle get boundary points from
            % t = tan(theta/2), canonical form
            %
            %  Inputs:
            %    t: Tangent of half angle, Nx2 vector
            %
            %  Outputs:
            %    x: Surface point, 3xN vector

            x = nan(3,size(t,1));
            
            [cos_t1, sin_t1] = trig_subs(t(:,1));
            [cos_t2, sin_t2] = trig_subs(t(:,2));
            
            x(1,:) = obj.a(1) .* eps_fun(cos_t1, obj.eps(1))...
                .* eps_fun(cos_t2, obj.eps(2));
            x(2,:) = obj.a(2) .* eps_fun(cos_t1, obj.eps(1))...
                .* eps_fun(sin_t2, obj.eps(2));
            x(3,:) = obj.a(3) .* eps_fun(sin_t1, obj.eps(1));
        end
        
        %% Get surface gradient cloud
        function m_t = GetGradients(obj)
            % GetGradients Get the surface gradients in global frame
            %
            %  Outputs:
            %    m_t: Surface gradients in global frame, 3xN vector
            
            m = obj.GetGradientsCanonical();
            x = obj.GetPointsFromGradient(m);
            m_t = quat2rotm(obj.q) * obj.TaperGradient(m, x);
        end
        
        function m = GetGradientsCanonical(obj)
            % GetGradientsCanonical Get gradient vectors (un-normalized) on
            % boundary, in canonical form
            %
            %  Outputs:
            %    m: Surface gradients in body frame, 3xN vector
            
            [m_x, m_y, m_z] = obj.GetSurfGradientsCanonical();
            
            % Reshape to 3xN array
            [mm, nn] = size(m_x);
            m = nan(3,mm*nn);
            m(1,:) = reshape(m_x, 1, mm*nn);
            m(2,:) = reshape(m_y, 1, mm*nn);
            m(3,:) = reshape(m_z, 1, mm*nn);
        end
        
        function [m_x, m_y, m_z] = GetSurfGradientsCanonical(obj)
            % GetSurfGradientsCanonical Get gradient vectors
            % (un-normalized) on boundary, canonical form, as mesh grid
            %
            %  Outputs:
            %    m_x, m_y, m_z: Surface gradients, in the form of mesh grid
            
            m_x = 2/(obj.a(1) * obj.eps(1)) .*...
                eps_fun(cos(obj.eta), 2-obj.eps(1)) .*...
                eps_fun(cos(obj.omega), 2-obj.eps(2));
            m_y = 2/(obj.a(2) * obj.eps(1)) .*...
                eps_fun(cos(obj.eta), 2-obj.eps(1)) .*...
                eps_fun(sin(obj.omega), 2-obj.eps(2));
            m_z = 2/(obj.a(3) * obj.eps(1)) .*...
                eps_fun(sin(obj.eta), 2-obj.eps(1)) .*...
                ones(size(obj.omega));
        end
        
        %% Get gradients from different parameterizations
        function m = GetGradientsFromParamTapered(obj, par, opt)
            % GetGradientsFromParamTapered Get surface gradients from '
            % different parameterizations
            %
            %  Inputs:
            %    par: Parameters
            %    opt: Indicator of parameterization
            %
            %  Outputs:
            %    x  : Surface gradients, 3xN vector
            
            if strcmp(opt, 'spherical')
                m = obj.GetGradientsFromSpherical(par)';
                x = obj.GetPointsFromSpherical(par);
                
            elseif strcmp(opt, 'cartesian')
                m = obj.GetGradientsFromCartesian(par);
                x = par;
                
            elseif strcmp(opt, 'hypersphere')
                m = obj.GetGradientsFromHypersphere(par);
                x = obj.GetPointsFromGradient(m);
                
            elseif strcmp(opt, 'tan-half-angle')
                m = obj.GetGradientsFromTanHalfAngle(par);
                x = obj.GetPointsFromTanHalfAngle(par);

            elseif strcmp(opt, 'direction')
                m = obj.GetGradientsFromDirection(par);
                x = obj.GetPointsFromGradient(m);
                
            else
                error('No such parameterization for surface gradients!')
            end
            
            % Tapering deformation
            m = obj.TaperGradient(m, x);
        end
        
        function m = GetGradientsFromSpherical(obj, theta)
            % GetGradientsFromSpherical Get gradient (un-normalized) on
            % boundary from spherical parameters
            %
            %  Inputs:
            %    theta: Angular parameters, Nx2 vector
            %
            %  Outputs:
            %    m    : Gradients for angular params, 3xN vector
            
            x = obj.GetPointsFromSpherical(theta);
            m = obj.GetGradientsFromCartesian(x);
        end
            
        function m = GetGradientsFromCartesian(obj, x)
            % GetGradientsFromCartesian Get surface gadient from cartisian
            % coordinates
            % 
            %  Inputs:
            %    x: Cartesian coordinates of points, 3xN vector
            %
            %  Outputs:
            %    m: Gradients of the points, 3xN vector
            
            m = nan(3,size(x,2));
            
            g_x1x2 = eps_fun( (x(1,:)./obj.a(1)).^2, 1/obj.eps(2) ) +...
                eps_fun( (x(2,:)./obj.a(2)).^2, 1/obj.eps(2) );
            
            m(1,:) = 2/(obj.a(1)*obj.eps(1)) *...
                eps_fun( g_x1x2, obj.eps(2)/obj.eps(1)-1 ) .*...
                eps_fun(x(1,:)/obj.a(1), 2/obj.eps(2)-1);
            m(2,:) = 2/(obj.a(2)*obj.eps(1)) *...
                eps_fun( g_x1x2, obj.eps(2)/obj.eps(1)-1 ) .*...
                eps_fun(x(2,:)/obj.a(2), 2/obj.eps(2)-1);
            m(3,:) = 2/(obj.a(3)*obj.eps(1)) *...
                eps_fun(x(3,:)/obj.a(3), 2/obj.eps(1)-1);
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
            
            m = nan(3,size(u,2));
            
            m(1,:) = 2/(obj.a(1)*obj.eps(1)) .*...
                eps_fun( u(1,:), 2-obj.eps(2) ) .*...
                eps_fun( u(1,:).^2+u(2,:).^2, (obj.eps(2)-obj.eps(1))/2 );
            m(2,:) = 2/(obj.a(2)*obj.eps(1)) .*...
                eps_fun( u(2,:), 2-obj.eps(2) ) .*...
                eps_fun( u(1,:).^2+u(2,:).^2, (obj.eps(2)-obj.eps(1))/2 );
            m(3,:) = 2/(obj.a(3)*obj.eps(1)) .*...
                eps_fun( u(3,:), 2-obj.eps(1) );
        end
        
        function m = GetGradientsFromTanHalfAngle(obj, t)
            % GetGradientsFromSpherical Get gradient (un-normalized) on
            % boundary from t = tan(theta/2)
            %
            %  Inputs:
            %    t: Tangent of half angle, Nx2 vector
            %
            %  Outputs:
            %    m: Gradients for angular params, 3xN vector
            
            x = obj.GetPointsFromTanHalfAngle(t);
            m = obj.GetGradientsFromCartesian(x);
        end
        
        function m = GetGradientsFromDirection(obj, v)
            % GetGradientsFromDirection Get gradient along any direction v
            %
            %  Inputs:
            %    v: Vectors indicating the direction, 3xN vector
            %
            %  Outputs:
            %    m: Gradients, 3xN vector
            
            dualObj = SuperQuadrics({2./(obj.a*obj.eps(1)),...
                2-obj.eps, [0,0], obj.tc, obj.q, obj.N});
            dualImplicit = dualObj.GetImplicitFunction(v);
            m = v./eps_fun(dualImplicit+1, dualObj.eps(1)/2);
        end
        
        %% Get surface normal cloud
        function n = GetNormals(obj)
            % GetNormals Get normal vectors on boundary
            %
            %  Outputs:
            %    n: Normal vectors, 3xN vector
            
            m = obj.GetGradients();
            
            % Normalize
            n = m./sqrt(sum(m.^2,1));
        end
        
        function [n_x, n_y, n_z] = GetNormalsCanonical(obj)
            % GetNormalsCanonical Get normal vectors on boundary, canonical
            % form, as mesh grid
            %
            %  Outputs:
            %    n_x, n_y, n_z: Surface normals, in the form of mesh grid
            
            [m_x, m_y, m_z] = obj.GetGradientsCanonical();
            
            gradNorm = sqrt( m_x.^2 + m_y.^2 + m_z.^2 );
            
            n_x = m_x ./ gradNorm;
            n_y = m_y ./ gradNorm;
            n_z = m_z ./ gradNorm;
        end
        
        function n_t = GetNormalsTransformed(obj)
            % GetNormalsTransformed Get normal at boundary points in global
            % frame
            %
            %  Outputs:
            %    n_t: Surface normals in global frame, 3xN vector
            
            n_t = quat2rotm(obj.q) * obj.GetNormals();
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
            %    u: Hypersphere parameters, 3xN vector
            
            u = nan(3,size(m,2));
            
            u(1,:) = eps_fun( (obj.a(1)*obj.eps(1))/2 .* m(1,:),...
                1/(2-obj.eps(2)) ) .*...
                eps_fun( eps_fun( (obj.a(1)*obj.eps(1)/2*m(1,:)).^2,...
                1/(2-obj.eps(2))) +...
                eps_fun( (obj.a(2)*obj.eps(1)/2*m(2,:)).^2,...
                1/(2-obj.eps(2))),...
                (obj.eps(1)-obj.eps(2))/(4-2*obj.eps(1)) );
            u(2,:) = eps_fun( (obj.a(2)*obj.eps(1))/2 .* m(2,:),...
                1/(2-obj.eps(2)) ) .*...
                eps_fun( eps_fun( (obj.a(1)*obj.eps(1)/2*m(1,:)).^2,...
                1/(2-obj.eps(2))) +...
                eps_fun( (obj.a(2)*obj.eps(1)/2*m(2,:)).^2,...
                1/(2-obj.eps(2))),...
                (obj.eps(1)-obj.eps(2))/(4-2*obj.eps(1)) );
            u(3,:) = eps_fun( (obj.a(3)*obj.eps(1))/2 .* m(3,:),...
                1/(2-obj.eps(1)) );
        end
        
        %% Deformations
        function x = TaperPoint(obj, x)
            % TaperPoint Get boundary points after tapering deformation
            %
            %  Inputs:
            %    x: Original boundary points, 3xN vector
            %
            %  Outputs:
            %    x: Tapered boundary points, 3xN vector
            
            x(1,:) = (obj.taper(1) / obj.a(3) .* x(3,:) + 1) .* x(1,:);
            x(2,:) = (obj.taper(2) / obj.a(3) .* x(3,:) + 1) .* x(2,:);
        end
        
        function [x, y, z] = TaperSurfPoints(obj, x, y, z)
            % TaperSurfPoints Get boundary points after tapering 
            % deformation, as mesh grid
            %
            %  Inputs:
            %    x, y, z: Original boundary points, as mesh grid
            %
            %  Outputs:
            %    x, y, z: Tapered boundary points, as mesh grid
            
            x = (obj.taper(1) / obj.a(3) .* z + 1) .* x;
            y = (obj.taper(2) / obj.a(3) .* z + 1) .* y;
        end
        
        function x = DeTaperPoint(obj, x)
            % DeTaperPoint Get original boundary points from tapering 
            % deformation
            %
            %  Inputs:
            %    x: Tapered boundary points, 3xN vector
            %
            %  Outputs:
            %    x: Original boundary points, 3xN vector
            
            x(1,:) = x(1,:) ./ (obj.taper(1)/obj.a(3) .* x(3,:) + 1);
            x(2,:) = x(2,:) ./ (obj.taper(2)/obj.a(3) .* x(3,:) + 1);
        end
        
        function m = TaperGradient(obj, m, x)
            % TaperGradient Get boundary point gradient after tapering
            % deformation
            %
            %  Inputs:
            %    m: Original gradients of boundary points, 3xN vector
            %    x: Original boundary points, 3xN vector
            %
            %  Outputs:
            %    m: Tapered gradients, 3xN vector
            
            fx = (obj.taper(1) / obj.a(3) .* x(3,:) + 1);
            fy = (obj.taper(2) / obj.a(3) .* x(3,:) + 1);
            fxp = obj.taper(1) / obj.a(3);
            fyp = obj.taper(2) / obj.a(3);
            
            m(1,:) = fy .* m(1,:);
            m(2,:) = fx .* m(2,:);
            m(3,:) = -fxp .* fy .* x(1,:) .* m(1,:) -...
                fyp .* fx .* x(2,:) .* m(2,:) +...
                fx .* fy .* m(3,:);
        end
        
        %% GetSurf Get mesh information of superquadrics, as surface object
        function S = GetSurf(obj)
            % GetSurf Get surface object of the superquadrics
            %
            %  Outputs:
            %    S: a "surf" object
            %
            %  See also
            %    surf
            
            [x, y, z] = obj.GetSurfPoints();
            
            S = surf(x, y, z, 'FaceColor', 'none', 'EdgeColor', 'none');
        end
        
        %% Plot superquadric surface
        function PlotShape(obj, color, faceAlpha)
            % PlotShape Plot the superquadric surface
            %
            %  Inputs:
            %    color    : color of the surface
            %    faceAlpha: transparency of the surface, range in [0,1]
            
            [x, y, z] = obj.GetSurfPoints();
            
            surf(x, y, z, 'FaceColor', color, 'EdgeColor', 'none',...
                'FaceAlpha', faceAlpha);
        end
        
        %% GetVolume compute volume of superquadrics
        function vol = GetVolume(obj)
            % GetVolume Compute the volume of the superquadrics
            %
            %  Outputs:
            %    vol: volume of the superquadric model
            
            vol = 2 * prod(obj.a) * prod(obj.eps) *...
                beta(obj.eps(1)/2+1, obj.eps(1)) *...
                beta(obj.eps(2)/2, obj.eps(2)/2);
        end
    end
end