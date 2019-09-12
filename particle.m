
% Particle member functions for correspondence likelihood.
%
classdef particle < handle
    
    properties
        pose = [0.0,0.0,0.0];
        landmark_positions 
        landmark_covariances
        counter = 0;
    end %properties
  
    methods         
      function self = particle(pose,landmark_positions,landmark_covariances,counter)
          if nargin > 0
              if(isa(pose,'particle'))
                  self.pose = pose.pose;
                  self.landmark_positions(:,:,1) = pose.landmark_positions;
                  self.landmark_covariances(:,:,1) = pose.landmark_covariances;
                  self.counter = pose.counter;
              else
                  self.pose = pose;
                  self.landmark_positions(:,:,1) = landmark_positions;
                  self.landmark_covariances(:,:,1) = landmark_covariances;
                  self.counter = counter;
              end
          end
      end %function constructor particle
      
      function total_landmarks = number_of_landmarks(self)
        % The size of self is 1 x 1
        % return current number of landmarks in this particle.%
        total_landmarks = length(self(1).landmark_positions);
      end %function number_of_landmarks

      function self = move(self, left, right, w)
        % The size of self is 1 x 1
        % Given left, right control and robot width, move the robot.
        state = self(1).pose;
        self(1).pose = self(1).stateTransition(state, [left, right], w);
      end %function move

      function [range, bearing] = h(self, landmark_number, scanner_displacement)
        % The size of self is 1 x 1
        %Measurement function. Takes a (x, y, theta) state and a (x, y)
        %landmark, and returns the corresponding (range, bearing).%
        state = self.pose;
        landmark = self.landmark_positions(:,:,landmark_number);
        dx = landmark(1) - (state(1) + scanner_displacement * cos(state(3)));
        dy = landmark(2) - (state(2) + scanner_displacement * sin(state(3)));
        range = sqrt(dx * dx + dy * dy);
        bearing = mod((atan2(dy, dx) - state(3) + pi),(2*pi)) - pi;
      end %function h
      
      function [dh_dlandmark] = dh_dlandmark(self, landmark_number, scanner_displacement)
        % The size of self is 1 x 1
        %Derivative with respect to the landmark coordinates.
        state = self.pose;
%         landmark_number is used as a switch as it can be used as landmark
%         container or as an index of landmark position
        if (length(landmark_number)==1)
            landmark = self.landmark_positions(:,:,landmark_number);
        elseif(length(landmark_number)==2)
            landmark = landmark_number;
        else
            display(' ERROR Length of Variable Landmakr Number is neighter 1 nor 2')
            landmark = [0,0];
        end
        theta = state(3);
        cost = cos(theta); 
        sint = sin(theta);
        dx = landmark(1) - (state(1) + scanner_displacement * cost);
        dy = landmark(2) - (state(2) + scanner_displacement * sint);
        q = dx * dx + dy * dy;
        sqrtq = sqrt(q);
        dr_dmx = dx / sqrtq;
        dr_dmy = dy / sqrtq;
        dalpha_dmx = -dy / q;
        dalpha_dmy =  dx / q;
        dh_dlandmark = [dr_dmx, dr_dmy;dalpha_dmx, dalpha_dmy];
      end %function dh_dlandmark
  
      function [expected_distance, expected_angle] = h_expected_measurement_for_landmark(self, landmark_number, scanner_displacement)
        % The size of self is 1 x 1
        % Returns the expected distance and bearing measurement for a given Landmark number and the pose of this particle.%
        % - the state is the robot's pose.
        state = self.pose;
        % - the landmark is taken from self.landmark_positions.
        landmark = self.landmark_positions(:,:,landmark_number);
        % - the static function h() computes the desired value
        [expected_distance, expected_angle] = h(self, landmark_number, scanner_displacement);

      end %function h_expected_measurement_for_landmark
      
      function [H,Ql] = H_Ql_jacobian_and_measurement_covariance_for_landmark(self, landmark_number, Qt_measurement_covariance, scanner_displacement)
        % The size of self is 1 x 1
        %Computes Jacobian H of measurement function at the particle's position and the landmark given by landmark_number. 
        %Also computes the measurement covariance matrix.%
        % - H is computed using dh_dlandmark.
        H = dh_dlandmark(self,landmark_number, scanner_displacement);
        landmark_cov = self.landmark_covariances(:,:,landmark_number);
        Ql = H * landmark_cov * H' + Qt_measurement_covariance;
      end %function H_Ql_jacobian_and_measurement_covariance_for_landmark
      
      function measurement_landmark_likelihood_correspondence = wl_likelihood_of_correspondence(self, measurement, landmark_number, Qt_measurement_covariance, scanner_displacement)
        % The size of self is 1 x 1
        %For a given measurement and landmark_number, returns the likelihood
        % 1 x 1 The size of self in h_expected_measurement_for_landmark()
        [expected_distance, expected_angle] = h_expected_measurement_for_landmark(self, landmark_number,scanner_displacement);
        expected_measurement = [expected_distance, expected_angle];
        % deltaz must be 1x2
        deltaZ = [measurement - expected_measurement]'; 
        %   will only need Ql, not H
        % 1 x 1 The size of self in H_Ql_jacobian_and_measurement_covariance_for_landmark
        [H,Ql] = H_Ql_jacobian_and_measurement_covariance_for_landmark(self, landmark_number, Qt_measurement_covariance, scanner_displacement);
        Ql_inv = inv(Ql);
        val = -.5 * deltaZ'* Ql_inv *deltaZ;
        detQl = det(Ql);
        measurement_landmark_likelihood_correspondence = 1 / (2 * pi * sqrt(detQl)) * exp(val);
      end %function wl_likelihood_of_correspondence
      
      function likelihoods = compute_correspondence_likelihoods(self, measurement, number_of_landmarks, Qt_measurement_covariance, scanner_displacement)
        % The size of self is 1 x 1
        % For a given measurement, returns a list of all correspondence likelihoods (from index 1 to number_of_landmarks).%
        likelihoods = [];
        if number_of_landmarks > 0
            for i = 1:number_of_landmarks
                % 1 x 1 The size of self in
                % wl_likelihood_of_correspondence()
                likelihoods = [likelihoods; wl_likelihood_of_correspondence(self,measurement, i, Qt_measurement_covariance, scanner_displacement)];
            end %for i = 1:number_of_landmarks
        end %if number_of_landmarks > 0
      end %function compute_correspondence_likelihoods
      
      function self = initialize_new_landmark(self, measurement_in_scanner_system,Qt_measurement_covariance,scanner_displacement)
        % The size of self is 1 x 1
        %Given a (x, y) measurement in the scanner's system, initializes a new landmark and its covariance.%
         state = self.pose;
         scanner_pose = [state(1) + cos(state(3)) * scanner_displacement,...
                        state(2) + sin(state(3)) * scanner_displacement,...
                        state(3)];
        [mx,my] = self(1).scanner_to_world(scanner_pose, measurement_in_scanner_system);
        landmark_positions = [mx,my];
         % - H is obtained from dh_dlandmark()
        dh_dlandmark1 = dh_dlandmark(self, landmark_positions, scanner_displacement);
        % - Use inv(A) to invert matrix A
        dh_dlandmark_inv = inv(dh_dlandmark1);
        landmark_covariances = dh_dlandmark_inv * Qt_measurement_covariance * dh_dlandmark_inv';
        
        if(isempty(self.landmark_positions))
            self.landmark_positions(:,:,end) = landmark_positions;
            self.landmark_covariances(:,:,end)=landmark_covariances;
            self.counter = self.counter + 1;
        else
            self.landmark_positions(:,:,end+1) = landmark_positions;
            self.landmark_covariances(:,:,end+1)=landmark_covariances;
            self.counter = self.counter + 1;   
        end %if(isempty(self.landmark_positions))
      end %function  initialize_new_landmark
      
      function  self = update_landmark(self, landmark_number, measurement, Qt_measurement_covariance, scanner_displacement)
        % The size of self is 1 x 1
        % Update a landmark's estimated position and covariance.%
        old_landmark_positions = self.landmark_positions(:,:,landmark_number);
        old_landmark_covariances = self.landmark_covariances(:,:,landmark_number);
        [H,Ql] = H_Ql_jacobian_and_measurement_covariance_for_landmark(self, landmark_number, Qt_measurement_covariance, scanner_displacement);
        Ql_inv = inv(Ql);
        % Compute Kalman gain
        K = old_landmark_covariances * H' * Ql_inv;
        % - Expected measurement can be computed using h_expected_measurement_for_landmark()
        [expected_distance, expected_angle] = h_expected_measurement_for_landmark(self, landmark_number, scanner_displacement);
        % - Delta z is measurement minus expected measurement
        delta_z = measurement-[expected_distance,expected_angle];
        % - update landmark_positions[landmark_number]
        new_landmark_positions = old_landmark_positions + (K * delta_z')';
        % - update landmark_covariances[landmark_number].
        new_landmark_covariances = old_landmark_covariances - (K * H * old_landmark_covariances);
        self.landmark_positions(:,:,landmark_number) = new_landmark_positions;
        self.landmark_covariances(:,:,landmark_number) = new_landmark_covariances;
      end %function  update_landmark
      
      function [self,maximum_likelihood] = update_particle(self, measurement, measurement_in_scanner_system,number_of_landmarks,minimum_correspondence_likelihood,Qt_measurement_covariance, scanner_displacement)
        % The size of self is 1 x 1
        % Given a measurement, computes the likelihood that it belongs to any of the landmarks in the particle. If there are none, or if all
        % likelihoods are below the minimum_correspondence_likelihood threshold, add a landmark to the particle. Otherwise, update the
        % (existing) landmark with the largest likelihood. Compute likelihood of correspondence of measurement to all landmarks
        % (from 1 to number_of_landmarks).
        likelihoods = compute_correspondence_likelihoods(self, measurement, number_of_landmarks, Qt_measurement_covariance, scanner_displacement);
        % compute_correspondence_likelihoods().
        % If the likelihood list is empty, or the max correspondence likelihood is still smaller than minimum_correspondence_likelihood, setup a new landmark.
        if ((isempty(likelihoods)) || (max(likelihoods) < minimum_correspondence_likelihood))
            self= initialize_new_landmark(self, measurement_in_scanner_system,Qt_measurement_covariance,scanner_displacement);
            maximum_likelihood = minimum_correspondence_likelihood;
        % Else update the particle's EKF for the corresponding particle.
        else
            % This computes (max, argmax) of measurement_likelihoods.
            % -find w, the maximum_likelihood, and the corresponding landmark index.
            [maximum_likelihood, index] = max(likelihoods); 
            % Add code to update_landmark().
            self = update_landmark(self, index, measurement, Qt_measurement_covariance, scanner_displacement);
        end %if ((isempty(likelihoods)) || (max(likelihoods) < minimum_correspondence_likelihood))
      end %function update_particle
    end  %Methods
    
    methods (Static)
      function pose = stateTransition(pose, control, w)
        %State transition.
        x = pose(1); 
        y = pose(2); 
        theta = pose(3);
        l = control(1);  
        r = control(2);
        
        if (r ~= l)
            alpha = (r - l) / w;
            rad = l/alpha;
            g1 = x + (rad + w/2.)*(sin(theta+alpha) - sin(theta));
            g2 = y + (rad + w/2.)*(-cos(theta+alpha) + cos(theta));
            g3 = mod((theta + alpha + pi),(2*pi)) - pi;
            pose =[g1,g2,g3];
        else
            g1 = x + l * cos(theta);
            g2 = y + l * sin(theta);
            g3 = theta;
            pose =[g1,g2,g3];
        end % if (r ~= l) else
      end %function stateTransition
      function [worldx,worldy] = scanner_to_world(pose, point)
%         Given a robot pose (rx, ry, heading) and a point (x, y) in the
%         scanner's coordinate system, return the point's coordinates in the
%         world coordinate system.
        dx = cos(pose(3));
        dy = sin(pose(3));
        x = point(1);
        y = point(2);
        worldx = (x * dx - y * dy + pose(1)); 
        worldy = (x * dy + y * dx + pose(2));     
    end
    end %Methods Static
    
end %classdef particle



