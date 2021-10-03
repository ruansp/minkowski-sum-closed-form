classdef SuperEllipse < handle
    % SuperEllipse Class of Superellipse operations, a SuperQuadric in 2D
    % 
    %  Inputs:
    %    val(1:2): Semi-axis lengths
    %    val(3)  : Exponent power
    %    val(4)  : Tapering factor, default = 0
    %    val(5:6): Center of the superellipse
    %    val(7)  : Angle of the superellipse
    %    val(8)  : Number of points to be sampled on boundary
    %
    %  Author:
    %    Sipu Ruan, ruansp@nus.edu.sg, 2021
    
    properties
        a       % length of semi-axes
        eps     % exponent for the signed power function
        taper   % tapering factor
        
        % pose of the object
        tc      % center of superellipse
        ang     % angle for the orientation
        
        N       % number of vertices on the boundary
    end
    
    properties (Access = private)
        theta   % Angle parameters
    end
    
    methods
        %% Constructor
        function obj = SuperEllipse(val)
            if length(val) ~= 8
                error('Length of 1st input not equal to 8!')
            else
                obj.a     = val(1:2);
                obj.eps   = val(3);
                obj.taper = val(4);
                obj.tc    = val(5:6)';
                obj.ang   = val(7);
                
                obj.N = val(8);
                obj.theta = GenerateParameters(obj, obj.N);
            end
        end
        
        %% Set the number of vertices on boundary
        function obj = SetNumVertices(obj, N)
            % SetNumVertices Set the number of vertices to be sampled on
            % boundary
            %
            %  Inputs:
            %    N: Number of vertices to be sampled on boundary
            
            obj.theta = GenerateParameters(obj, N);
            obj.n_vtx = N;
        end
        
        %% Get implicit function value
        function f = GetImplicitFunction(obj, x)
            % GetImplicitFunction Get implicit expression value of
            % superellipse for a point in space. The point is expressed in
            % global frame
            %
            %  Inputs:
            %    x: A set of points in space, 2xN vector
            %
            %  Outputs:
            %    f: implicit expression value of superellipse, 1xN vector

            x = obj.DeTaperPoint(x);
            f = obj.GetImplicitFunctionCanonical(x);
        end
        
        function f = GetImplicitFunctionCanonical(obj, x)
            % GetImplicitFunctionCanonical Get implicit expression value of
            % superellipse for a point in space
            %
            %  Inputs:
            %    x: A set of points expressed in body frame, 2xN vector
            %
            %  Outputs:
            %    f: implicit expression value of superellipse, 1xN vector
            
            f = eps_fun( (x(1,:) ./ obj.a(1)).^2, 1/obj.eps ) +...
                eps_fun( (x(2,:) ./ obj.a(2)).^2, 1/obj.eps ) - 1;
        end
        
        %% Get points on boundary from different parameterizations
        function pnt = GetPoints(obj)
            % GetPoints Get points on boundary in global frame
            %
            %  Outputs:
            %    pnt: The point on boundary in global frame, 2xN vector
            
            pnt_canonical = obj.GetPointsCanonical();
            pnt_tapered = obj.TaperPoint(pnt_canonical);
            pnt = rot2(obj.ang) * pnt_tapered + obj.tc;
        end
        
        function pnt = GetPointsCanonical(obj)
            % GetPointsCanonical Get points on boundary in body frame
            %
            %  Outputs:
            %    pnt: The point on boundary in body frame, 2xN vector
            
            pnt(1,:) = obj.a(1) * eps_fun(cos(obj.theta), obj.eps);
            pnt(2,:) = obj.a(2) * eps_fun(sin(obj.theta), obj.eps);
        end
        
        function f_m = GetPointsFromParamTapered(obj, par, opt)
            % GetPointsFromParamTapered Get boundary points using
            % different parameterizations after global deformations
            %
            %  Inputs:
            %    par: Parameters
            %    opt: Indicator of parameterization
            %
            %  Outputs:
            %    f_m: A set of boundary points, 2xN vector
            
            if strcmp(opt, 'gradient')
                f_m = obj.GetPointsFromGradient(par);
            end
            
            % Tapering deformation
            f_m = obj.TaperPoint(f_m);
        end
        
        function f_m = GetPointsFromGradient(obj, m)
            % GetPointsFromGradient Get boundary points from gradient
            % vectors
            %
            %  Inputs:
            %    m  : A set of gradient vectors, 2xN vector
            %
            %  Outputs:
            %    f_m: cartesian coordinates of points, 2xN vector
            
            f_m(1,:) = obj.a(1) .* eps_fun( obj.eps*obj.a(1)/2 * m(1,:),...
                obj.eps/(2-obj.eps) );
            f_m(2,:) = obj.a(2) .* eps_fun( obj.eps*obj.a(2)/2 * m(2,:),...
                obj.eps/(2-obj.eps) );
        end
        
        %% Get gradient on boundary from different parameterizations
        function m_t = GetGradients(obj)
            % GetGradients Get gradient vectors in global frame
            %
            %  Outputs:
            %    m_t: A set of gradient vectors in global frame, 2xN vector
            
            m = obj.GetGradientsCanonical();
            x = obj.GetPointsFromGradient(m);
            m_t = rot2(obj.ang) * obj.TaperGradient(m, x);
        end
        
        function m = GetGradientsCanonical(obj)
            % GetGradientsCanonical Get gradient vectors in body frame
            %
            %  Outputs:
            %    m: A set of gradient vectors in body frame, 2xN vector
            
            m(1,:) = 2/(obj.a(1)*obj.eps) * ...
                eps_fun(cos(obj.theta), 2-obj.eps);
            m(2,:) = 2/(obj.a(2)*obj.eps) * ...
                eps_fun(sin(obj.theta), 2-obj.eps);
        end
        
        function m = GetGradientsFromParamTapered(obj, par, opt)
            % GetGradientsFromParamTapered Get gradient vector using
            % different parameterizations after global deformations
            %
            %  Inputs:
            %    par: Parameters
            %    opt: Indicator of parameterization
            %
            %  Outputs:
            %    m  : A set of gradient vectors, 2xN vector
            
            if strcmp(opt, 'angle')
                m = obj.GetGradientFromAngle(par);
                
            elseif strcmp(opt, 'cartesian')
                m = obj.GetGradientsFromCartesian(par);
                
            elseif strcmp(opt, 'hypersphere')
                m = obj.GetGradientsFromHypersphere(par);
            end
            
            % Tapering deformation
            x = obj.GetPointsFromGradient(m);
            m = obj.TaperGradient(m, x);
        end
        
        function m = GetGradientFromAngle(obj, theta)
            % GetGradientFromAngle Get gradient vectors from angular
            % parameters
            %
            %  Inputs:
            %    theta: A set of angular parameters, 1xN vector
            %
            %  Outputs:
            %    m    : A set of gradient vectors, 2xN vector
            
            m(1,:) = 2/(obj.a(1)*obj.eps) * eps_fun(cos(theta), 2-obj.eps);
            m(2,:) = 2/(obj.a(2)*obj.eps) * eps_fun(sin(theta), 2-obj.eps);
        end
        
        function m = GetGradientsFromCartesian(obj, x)
            % GetGradientsFromCartesian Get gradient vectors from boundary
            % points
            %
            %  Inputs:
            %    x: Cartesian coordinates of points, 2xN vector
            %
            %  Outputs:
            %    m: A set of gradient vectors, 2xN vector
            
            m(1,:) = 2/(obj.a(1)*obj.eps) *...
                eps_fun( x(1,:)/obj.a(1), 2/obj.eps-1 );
            m(2,:) = 2/(obj.a(2)*obj.eps) *...
                eps_fun( x(2,:)/obj.a(2), 2/obj.eps-1 );
        end
        
        function m = GetGradientsFromHypersphere(obj, u)
            % GetGradientsFromHypersphere Get gradient vectors from
            % hypersphere parameters
            %
            %  Inputs:
            %     u: Unit hypersphere parameters, 2xN vector
            %
            %  Outputs:
            %    m: A set of gradient vectors, 2xN vector
            
            m(1,:) = 2/(obj.a(1)*obj.eps) .* eps_fun(u(1,:), 2-obj.eps);
            m(2,:) = 2/(obj.a(2)*obj.eps) .* eps_fun(u(2,:), 2-obj.eps);
        end
        
        %% Get unit hypersphere parameters from gradient
        function u = GetHypersphereFromGradient(obj, m)
            % GetHypersphereFromGradient Get hypersphere parameters from
            % gradient vectors
            %
            %  Inputs:
            %    m: A set of gradient vectors, 2xN vector
            %
            %  Outputs:
            %    u: A set of unit hypersphere parameters, 2xN vector
            
            u(1,:) = eps_fun( (obj.a(1)*obj.eps)/2 .* m(1,:),...
                1/(2-obj.eps) );
            u(2,:) = eps_fun( (obj.a(2)*obj.eps)/2 .* m(2,:),...
                1/(2-obj.eps) );
        end
        
        %% Global tapering deformations
        function x = TaperPoint(obj, x)
            % TaperPoint Get boundary points after tapering deformation
            %
            %  Inputs:
            %    x: Original boundary points, 2xN vector
            %
            %  Outputs:
            %    x: Tapered boundary points, 2xN vector
            
            x(1,:) = (obj.taper / obj.a(2) .* x(2,:) + 1) .* x(1,:);
        end
        
        function x = DeTaperPoint(obj, x)
            % DeTaperPoint Get original boundary points from tapering 
            % deformation
            %
            %  Inputs:
            %    x: Tapered boundary points, 2xN vector
            %
            %  Outputs:
            %    x: Original boundary points, 2xN vector
            
            x(1,:) = x(1,:) ./ (obj.taper / obj.a(2) .* x(2,:) + 1);
        end
        
        function m = TaperGradient(obj, m, x)
            % TaperGradient Get boundary point gradient after tapering
            % deformation
            %
            %  Inputs:
            %    m: Original gradients of boundary points, 2xN vector
            %    x: Original boundary points, 2xN vector
            %
            %  Outputs:
            %    m: Tapered gradients, 2xN vector
            
            fx = obj.taper / obj.a(2) .* x(2,:) + 1;
            fxp = obj.taper / obj.a(2);
            
            m(2,:) = -fxp .* x(1,:) .* m(1,:) + fx .* m(2,:);
        end
        
        %% Plot superelliptical shape
        function PlotShape(obj, color, faceAlpha)
            % PlotShape Plot the superellipse shape
            %
            %  Inputs:
            %    color    : the color of the surface patch
            %    faceAlpha: transparency of the patch, range in [0,1]
            %
            %  See also
            %    patch
            
            if nargin == 1
                color = 'b';
            end
            
            if nargin == 2
                faceAlpha = 0.6;
            end
            
            xy = obj.GetPoints();
            patch(xy(1,:), xy(2,:), color, 'FaceAlpha', faceAlpha);
        end
        
        %% Generate almost equally spaced parameters
        function theta = GenerateParameters(obj, N)
            % GenerateParameters Generate angular parameters of the
            % boundary
            %
            %  Inputs:
            %    N    : Number of points to be sampled in boundary
            %
            %  Outputs:
            %    theta: Angular parameters
            
            theta_0 = linspace(-pi, pi, N);
            
            % Almost equally spaced samples
            theta = atan( sign(tan(theta_0)) .* ...
                abs(tan(theta_0)).^(1/obj.eps) );
            
            % Adjust range of omega from [-pi/2,pi/2] to [-pi,pi]
            theta(theta_0 > pi/2) = theta(theta_0 > pi/2) + pi;
            theta(theta_0 < -pi/2) = theta(theta_0 < -pi/2) - pi;
        end
    end
end