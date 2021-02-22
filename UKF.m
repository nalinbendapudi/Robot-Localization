classdef UKF < handle
    properties
        mu;             % (3x1)Pose Mean 
        Sigma;          % (3x3) Pose Covariance 
        gfun;           % (3x1)Motion Model Function 
        hfun;           % (2x1)Measruement Model Function
        M;              % (3x3)Motion model noise(dynamical and function of input)
        Q;              % (2x2)Sensor Noise
        kappa_g;        %scalar
        mu_pred;        %(3x1)
        Sigma_pred;     %(3x3)
        n;              %scalar
        X;
        w;
        Y;
        X_star_x;
    end
    
    methods
        function obj = UKF(sys, init)
            % motion model
            obj.gfun = sys.gfun;
            % measurement model
            obj.hfun = sys.hfun;
            % motion noise covariance
            obj.M = sys.M;
            % measurement noise covariance
            obj.Q = sys.Q;
            obj.kappa_g = init.kappa_g;
            % initial mean and covariance
            obj.mu = init.mu;
            obj.Sigma = init.Sigma;
        end
        
        function prediction(obj, u)
            %Augment
            mu_aug = [obj.mu', [0,0,0],[0,0]]'; %(8x1)
            Sigma_aug = blkdiag(obj.Sigma, obj.M(u), obj.Q);
            %Create Sigma Points
            sigma_point(obj, mu_aug, Sigma_aug, obj.kappa_g);   
            
            %Propagate Sigma Points through Non-linear map
            obj.X_star_x = zeros(numel(obj.mu),numel(obj.w));
            for i=1:numel(obj.w)
                obj.X_star_x(:,i) = obj.gfun(obj.X(1:3,i), u+obj.X(4:6,i));
            end
            
            %Predict Mean
            obj.mu_pred = sum(obj.w' .* obj.X_star_x, 2);
            obj.mu_pred(3) = wrapToPi(obj.mu_pred(3));
            %Predict Covariance
            obj.Sigma_pred = zeros(size(obj.mu,1));     
            for i=1:numel(obj.w)
                X_x_diff = obj.X_star_x(:,i)-obj.mu_pred;
                X_x_diff(3) = wrapToPi(X_x_diff(3));
                temp = obj.w(i) * (X_x_diff * X_x_diff');
                obj.Sigma_pred = obj.Sigma_pred + temp;
            end
        end
        
        function correction(obj, z)
            global FIELDINFO;
            landmark_x = FIELDINFO.MARKER_X_POS(z(3));
            landmark_y = FIELDINFO.MARKER_Y_POS(z(3));
            z_star=zeros(numel(z)-1,numel(obj.w));
            for i=1:numel(obj.w)
                z_star(:,i) = obj.hfun(landmark_x, landmark_y,obj.X_star_x(:,i)) + obj.X(7:8,i);
                z_star(1,i) =  wrapToPi(z_star(1,i));
            end
            %Calculate Innovation
            z_hat = sum(obj.w' .* z_star,2);
            z_hat(1) = wrapToPi(z_hat(1));
            S = zeros(size(z,1) - 1);
            Sigma_xz = zeros(3,2);
            for i=1:size(obj.w)
                diff1 = z_star(:,i) - z_hat;
                diff1(1) = wrapToPi(diff1(1));
                temp1 = obj.w(i) * (diff1*diff1');
                S = S + temp1;
                
                diff2 = obj.X_star_x(:,i) - obj.mu_pred;
                diff2(3) = wrapToPi(diff2(3));
                temp2 = obj.w(i) * (diff2*diff1');
                Sigma_xz = Sigma_xz + temp2;
            end
            %Calculate Kalman Gain
            K = Sigma_xz / S;
            z(1) = wrapToPi(z(1));
            v = z(1:2) - z_hat;
            obj.mu = obj.mu_pred + K*(v);
            obj.Sigma = obj.Sigma_pred - K*S*K';
        end
        
        function sigma_point(obj, mean, cov, kappa)
            obj.n = numel(mean);
            L = sqrt(obj.n + kappa) * chol(cov,'lower');
            obj.Y = mean(:,ones(1, numel(mean)));
            obj.X = [mean,obj.Y + L, obj.Y - L];
            obj.w = zeros(2 * obj.n + 1,1);
            obj.w(1) = kappa / (obj.n + kappa);
            obj.w(2:end) = 0.5 / (obj.n + kappa);
        end
    end
end