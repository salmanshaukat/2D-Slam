% This file contains helper functions for the SLAM.
classdef helperLib < FastSLAM
    
    methods 

        % Find the derivative in scan data, ignoring invalid measurements.
        function jumps = compute_derivative(helperLibObj,scan, min_dist)%
            jumps = [];
            for i =1:length(scan)-1%
                l = scan(i);
                r = scan(i+1);
                if ((l > min_dist) && (r > min_dist))%
                    derivative = (r - l) / 2.0;
                    jumps=[jumps; derivative];
                else%
                    jumps=[jumps; 0];
                end
            end
            jumps=[jumps; 0];
        end

        function cylinder_list = find_cylinders(helperLibObj,scan, scan_derivative, jump, min_dist)
        % For each area between a left falling edge and a right rising edge,
        % determine the average ray number and the average depth.
            cylinder_list = [];
            on_cylinder = false;
            sum_ray = 0;
            sum_depth = 0.0;
            rays = 0;

            for i = 1:length(scan_derivative)
                if (scan_derivative(i) < -jump)
                    % Start a new cylinder, independent of on_cylinder.
                    on_cylinder = true;
                    sum_ray = 0;
                    sum_depth = 0.0;
                    rays = 0;
                elseif (scan_derivative(i) > jump)
                    % Save cylinder if there was one.
                    if (on_cylinder && rays)
                        cylinder_list = [cylinder_list;[sum_ray/rays, sum_depth/rays]];
                    end
                    on_cylinder = false;
                % Always add point, if it is a valid measurement.
                elseif (scan(i) > min_dist)
                    sum_ray = sum_ray + i;
                    sum_depth = sum_depth + scan(i);
                    rays = rays + 1;
                end
            %return cylinder_list
            end
        end

        function angle = beam_index_to_angle(helperLibObj,beam_index)
            mounting_angle = -0.06981317007977318;
            %Convert a beam index to an angle, in radians."""
            angle = (beam_index - 330.0) * 0.006135923151543 + mounting_angle;
        end

        function [cylinders] = get_cylinders_from_scan(helperLibObj,scan, jump, min_dist, cylinder_offset)
        % Detects cylinders and computes range, bearing and cartesian coordinates
        % (in the scanner's coordinate system).
        % The result is modified from previous versions it returns a list of
        % (distance, bearing, x, y) in the scanner's coordinate system.
            der = compute_derivative(helperLibObj,scan, min_dist);
            cylinders = find_cylinders(helperLibObj,scan, der, jump, min_dist);
            distanceS = []; bearingS = []; xS = []; yS = [];
            j=size(cylinders);
            for i =1: j(1);
                % Compute the angle and distance measurements.
                bearing = beam_index_to_angle(helperLibObj,cylinders(i,1));
                distance = cylinders(i,2) + cylinder_offset;

                bearingS = [bearingS; bearing];
                distanceS = [distanceS; distance];
                % Compute x, y of cylinder in the scanner system.
                x = distance * cos(bearing);
                y = distance * sin(bearing);

                xS = [xS; x];
                yS = [yS; y];
        %     return result
            end
            cylinders=[distanceS, bearingS, xS, yS];
        end

        function meanParticle = get_mean(helperLibObj,particlesObj)
        %Compute mean position and heading from a given set of particles."""
        % Note this function would more likely be a part of FastSLAM or a base class
        % of FastSLAM. It has been moved here for the purpose of keeping the
        % FastSLAM class short in this tutorial.
            mean_x = 0.0;
            mean_y = 0.0;
            head_x = 0.0;
            head_y = 0.0;
            n = length(particlesObj); % total number of particals
            for k = 1: n
                p = particlesObj(k).pose;
                x = p(1);
                y = p(2);
                theta = p(3);
                mean_x = mean_x + x;
                mean_y = mean_y + y;
                head_x = head_x + cos(theta);
                head_y = head_y + sin(theta);
            end

            meanParticle =[mean_x / n, mean_y / n, atan2(head_y, head_x)];
        end

        function [ellipse_angle,eigenvals0,eigenvals1,var_heading]=get_error_ellipse_and_heading_variance(particlesObj, mean)

            %{
            Given a set of particles and their mean (computed by get_mean()),
            returns a tuple: (angle, stddev1, stddev2, heading-stddev) which is
            the orientation of the xy error ellipse, the half axis 1, half axis 2,
            and the standard deviation of the heading.
            %}

            center_x  = mean(1);
            center_y  = mean(2);
            center_heading = mean(3);
            n = length(particlesObj);

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

                for k = 1: length(particlesObj)
                    p = particlesObj(k).pose;
                    x = p.pose(1);
                    y = p.pose(2);
                    theta = p.pose(3);
                    dx = x - center_x;
                    dy = y - center_y;
                    sxx = sxx + dx * dx;
                    sxy = sxy + dx * dy;
                    syy = syy + dy * dy;
                end

                cov_xy = [sxx, sxy; sxy, syy] / (n-1);

                % Get variance of heading.
                var_heading = 0.0;

                for k = 1: length(particlesObj)
                    p = particlesObj(k).pose;
                    dh = (p.pose(2) - center_heading + pi); % (2*pi) - pi
                    var_heading = var_heading + (dh * dh);
                end

                var_heading = var_heading / (n-1);

                % Convert xy to error ellipse.
                [eigenvals, eigenvects] = eig(cov_xy);

                % the rotation angle if the points in eigenvects are regarded as having 
                % been rotated counterclockwise to new positions through the angle theta
                % theta = atan2(V(2,1),V(1,1));

                ellipse_angle = atan2(eigenvects(2,1), eigenvects(1,1));
                eigenvals0 = (abs(eigenvals(1)));
                eigenvals1 = sqrt(abs(eigenvals(2)));
                var_heading = sqrt(var_heading);

            end
        end
        
        function allParticle = getSingleValue4rmAll(helperLibObj,FastSLAMobj)
        %Returns a double array containing position and heading from a given
        %set of particles.
            n = length(FastSLAMobj); % total number of particals
            allParticle = zeros(n,3);
            for k = 1: n
                p = FastSLAMobj(k).pose;
                x = p(1);
                y = p(2);
                theta = p(3);
                allParticle(k,:) = [x,y,theta];
            end
        end
                
    end%methods
end%classdef
