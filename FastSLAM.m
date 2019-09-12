classdef FastSLAM < particle 

    properties
        % Some constants.
        robot_width
        scanner_displacement
        control_motion_factor
        control_turn_factor
        measurement_distance_stddev
        measurement_angle_stddev
        minimum_correspondence_likelihood
    end % properties

    methods
        
        function self = FastSLAM(robot_width, scanner_displacement,control_motion_factor, control_turn_factor,measurement_distance_stddev, measurement_angle_stddev,minimum_correspondence_likelihood,varargin)
            self=self@particle(varargin{:});
            if nargin > 0
            self.robot_width = robot_width;
            self.scanner_displacement = scanner_displacement;
            self.control_motion_factor = control_motion_factor;
            self.control_turn_factor = control_turn_factor;
            self.measurement_distance_stddev = measurement_distance_stddev;
            self.measurement_angle_stddev = measurement_angle_stddev;
            self.minimum_correspondence_likelihood = minimum_correspondence_likelihood;
            end
        end %function FastSLAM constructor

        function self = predict(self, control)
            % The size of self can be 1 x N
            %The prediction step of the particle filter.
            left = control(1);  
            right = control(2);
            left_std  = sqrt((self(1).control_motion_factor * left)^2 +(self(1).control_turn_factor * (left-right))^2);
            right_std = sqrt((self(1).control_motion_factor * right)^2 +(self(1).control_turn_factor * (left-right))^2);
            % Modify list of particles: for each particle, predict its new position.
            for k = 1: length(self)
                l = random('Normal',left, left_std);
                r = random('Normal',right, right_std);
                p = self(k);
                p = move(p,l, r, self(1).robot_width);
                self(k) = p;
            end %for k = 1: length(self)
            
        end %function predict

        function [self, weights] = update_and_compute_weights(self, cylinders)
            %Updates all particles and returns a list of their weights.%
            % The size of self can be 1 x N
            memt_dist_stdev = self.measurement_distance_stddev;
            memt_angle_stdev = self.measurement_angle_stddev;
            min_corspnd_likelihood = self.minimum_correspondence_likelihood;
            scanr_dispcmt = self.scanner_displacement;
            Qt_measurement_covariance = diag([memt_dist_stdev^2, memt_angle_stdev^2]);
            weights = [];
            n = length(self); % total number of particals
            for k = 1: n
                number_of_landmarks =self(k).counter;
                weight = 1.0;
                if (~isempty(cylinders))
                    j=size(cylinders);
                    for i=1:j(1)
                        measurement = cylinders(i,1:2);
                        measurement_in_scanner_system =  cylinders(i,3:4);
                        % 1 x 1 The size of partical2Obj in update_particle
                        [self(k),maximum_likelihood] = update_particle(self(k),measurement, measurement_in_scanner_system,number_of_landmarks,min_corspnd_likelihood,Qt_measurement_covariance, scanr_dispcmt);
                        weight = weight * maximum_likelihood;
                    end %for i=1:j(1)
                end %if (~isempty(cylinders))
                % Append overall weight of this particle to weight list.
                weights = [weights, weight];
            end %for k = 1: n
        end %function update_and_compute_weights

        function self = resample(self, weights)
        % The size of partical2Obj is 1 x N
        %Return a list of particles which have been resampled, proportional %to the given weights.%
            new_particles_pose = []; new_particles_landmark_positions = {};  new_particles_landmark_covariances = {};  new_particles_counter = []; 
            max_weight = max(weights);
            index = randi([1,length(self)]);
            offset = 0.0;
            loop = length(self);
            for i =1:loop
                offset = offset + rand() * 2.0 * max_weight;
                weights(index);
                while (offset > weights(index))
                    offset =offset - weights(index);
                    index= index + 1;
                    index = round(mod((index), length(weights)));
                    if index > length(self)
                        index=length(self);
                    end %if index > length(self)
                    if index == 0
                        index = 1;
                    end %if index == 0
                end %while (offset > weights(index))
                if(i==1)
                    new_particles_pose(:,:,end) = self(index).pose;
                    new_particles_landmark_positions(i) = {self(index).landmark_positions};
                    new_particles_landmark_covariances(i) = {self(index).landmark_covariances};
                    new_particles_counter = self(index).counter;
                else
                    new_particles_pose(:,:,end+1) = self(index).pose;
                    new_particles_landmark_positions(i) = {self(index).landmark_positions};
                    new_particles_landmark_covariances(i) = {self(index).landmark_covariances};
                    new_particles_counter(:,end+1) = self(index).counter;
                end %if(i==1) else
            end %for i =1:loop
            for i =1:loop
                self(i).pose = new_particles_pose(:,:,i);
                self(i).landmark_positions = new_particles_landmark_positions{i}(:,:,:);
                self(i).landmark_covariances = new_particles_landmark_covariances{i}(:,:,:);
                self(i).counter = new_particles_counter(:,i);
            end %for i =1:loop
        end %function resample

        function self = correct(self, cylinders)
        %The correction step of FastSLAM.%
            % Update all particles and compute their weights.
            [self, weights] = update_and_compute_weights(self,cylinders);
            % Then resample, based on the weight array.
            % The size of self is 1 x N
            self= resample(self,weights);
        end %function correct
        
        function [ellipse_angle,eigenvals0,eigenvals1,var_heading]=get_error_ellipse_and_heading_variance(self, mean)
        % Given a set of particles and their mean (computed by get_mean()), returns a tuple: (angle, stddev1, stddev2, heading-stddev) which is
        %the orientation of the xy error ellipse, the half axis 1, half axis 2, and the standard deviation of the heading.
            center_x  = mean(1);
            center_y  = mean(2);
            center_heading = mean(3);
            n = length(self);
            if (n < 2)
                ellipse_angle = 0;
                eigenvals0 = 0;
                eigenvals1 = 0;
                var_heading = 0;
            else
            % Compute covariance matrix in xy.
            sxx = 0.0;
            sxy = 0.0;
            syy = 0.0;
            for k = 1: length(self)
                p = self(k).pose;
                x = p(1);
                y = p(2);
                theta = p(3);
                dx = x - center_x;
                dy = y - center_y;
                sxx = sxx + dx * dx;
                sxy = sxy + dx * dy;
                syy = syy + dy * dy;
            end
            cov_xy = [sxx, sxy; sxy, syy] / (n-1);
            % Get variance of heading.
            var_heading = 0.0;
            for k = 1: length(self)
                p = self(k).pose;
                dh = (p(2) - center_heading + pi); % (2*pi) - pi
                var_heading = var_heading + (dh * dh);
            end
            var_heading = var_heading / (n-1);
            % Convert xy to error ellipse.
            [eigenvals, eigenvects] = eig(cov_xy);
            % the rotation angle corresponds to the points if eigenvects are
            % assumed to have a counterclockwise rotation to new positions through the angle theta.
            % theta = atan2(V(2,1),V(1,1));
            ellipse_angle = atan2(eigenvects(2,1), eigenvects(1,1));
            eigenvals0 = (abs(eigenvals(1)));
            eigenvals1 = sqrt(abs(eigenvals(2)));
            var_heading = sqrt(var_heading);
            end
        end

    end
end