classdef InEKF < handle   
    properties
        mu;                 % Pose Mean
        Sigma;              % Pose Sigma
        gfun;               % Motion model function
        mu_pred;             % Mean after prediction step
        Sigma_pred;          % Sigma after prediction step
        mu_cart;            % Mean in cartesian coordinates
        sigma_cart;         % Sigma in cartesian coordinates
        N                   % Sensor Noise
        Q                   % Control Input Noise
    end
    
    methods
        function obj = InEKF(sys, init)
            obj.gfun = sys.gfun;
            obj.mu = init.mu;
            obj.Sigma = init.Sigma;
            obj.N = diag([1000, 1000]);            % sensor noise
            obj.Q = diag([0.015^2, 0.01^2, 0.01^2]); % control noise
        end
        
        function prediction(obj, u)
            
            % calculate mu in cartesian coordinates
            obj.mu_cart = [obj.mu(1,3); obj.mu(2,3); atan2(obj.mu(2,1),obj.mu(1,1))];
            
            % calculate twist corresponding to the motion model
            mu_pred_cart_temp = obj.gfun(obj.mu_cart, u);
            %mu_pred_cart_temp(3) = wrapTo2Pi(mu_pred_cart_temp(3));
            mu_pred_temp = obj.posemat(mu_pred_cart_temp);
            twist_hat = logm(obj.mu\mu_pred_temp);
            
            %Calculate adjoint of mu
            AdjX = [obj.mu(1:2,1:2), [obj.mu(2,3);-obj.mu(1,3)]; 0 0 1];
            
            propagation(obj, twist_hat, AdjX);
        
        end
        
        function propagation(obj, twist_hat, AdjX)
            
            % Adding noise to control input
            twist_noise = chol(obj.Q, 'lower')*rand(3,1);
            twist_noise_hat = wedge(obj, twist_noise);
            twist_hat = twist_hat + twist_noise_hat;
            
            % Prediction Update
            obj.mu_pred = obj.mu*expm(twist_hat);
            obj.Sigma_pred = obj.Sigma + AdjX*obj.Q*AdjX';
            
        end
        
        function correction(obj, Y1, Y2, landmark_ids)
            
            % get landmark coordinates
            global FIELDINFO;        
            landmark_x1 = FIELDINFO.MARKER_X_POS(landmark_ids(1));
            landmark_y1 = FIELDINFO.MARKER_Y_POS(landmark_ids(1));       
            landmark_x2 = FIELDINFO.MARKER_X_POS(landmark_ids(2));
            landmark_y2 = FIELDINFO.MARKER_Y_POS(landmark_ids(2));
            
            b1 = [landmark_x1;
                  landmark_y1;
                            1];
            b2 = [landmark_x2;
                  landmark_y2;
                            1];
            H1 = [-1 0  landmark_y1 ; 
                  0 -1 -landmark_x1 ];
            H2 = [-1 0  landmark_y2 ; 
                  0 -1 -landmark_x2 ];
            H = [H1;H2];
            
            % Adding sensor noise
            N_temp = obj.mu_pred*blkdiag(obj.N,0)*obj.mu_pred'; 
            N_temp = N_temp(1:2,1:2);
            S_Noise = blkdiag(N_temp,N_temp);
            
            % Calculating Kalman Gain
            S = H*obj.Sigma_pred*H' + S_Noise;
            L = (obj.Sigma_pred*H')*(S\eye(length(S))) ;
            
            % Calculating innovation
            nu1_cart = obj.mu_pred * Y1' - b1; 
            nu1_cart = nu1_cart(1:2);
            nu2_cart = obj.mu_pred * Y2' - b2; 
            nu2_cart = nu2_cart(1:2);
            nu_cart = [nu1_cart; nu2_cart];
            nu = obj.wedge(L*nu_cart);
            
            % Measurement Update
            obj.mu = expm(nu) * obj.mu_pred;
            I = eye(size(obj.Sigma_pred));
            obj.Sigma = (I-L*H)*obj.Sigma_pred*(I-L*H)' + L*S_Noise*L';            
            
        end
        
        function H = posemat(obj, state)
            x = state(1);
            y = state(2);
            h = state(3);
            % construct a SE(2) matrix element
            H = [...
                cos(h) -sin(h) x;
                sin(h)  cos(h) y;
                     0       0 1];
        end
        
        function xhat = wedge(obj,x)
            % converts R3 to se(2)
            % R3 should be (omega,vx,vy)
            xhat = [   0  -x(3) x(1);
                     x(3)    0  x(2);
                       0     0    0];
        end     
    end
end