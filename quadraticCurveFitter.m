clear
close all
clc

% Question 1
fprintf("Question 1\n\n\n") % First question

Yi = [202.36 239.03 280.71 309.12 323.15 332.78 328.45 306.40 287.36 247.97 202.89 161.11 93.68 20.78];
xi = 0:1:13;  % Dependent and independent values are entered

n = 2;% Highest degree in function is 2 since function is quadratic
A = zeros(n+1,n+1); % Sizes of the matrices are entered according to calculations on report
B = zeros(n+1,1);
x = zeros(n+1,1); % This part consists of coefficients of the function

for ii = n+1:-1:1 % Column number
    for jj = n+1:-1:1  % Row number
 % Row number increases for each column number until A and B matrices are
 % filled according to the report
        A(n+2-ii,n+2-jj)=sum(xi.^(ii+jj-2));
    end
    B(n+2-ii) = sum(xi.^(ii-1).*Yi);
end

x = A \ B;
            % Both equations give the same result
%x = inv(A)*B;

fprintf("The interpolated function is: f(x) = ") % Printing interpolated function
for ii=3:-1:2
    fprintf("%f x^%d + ",x(length(x)+1-ii),ii-1)    
end
fprintf("%f\n\n",x(length(x)))    


% Question 2
fprintf("Question 2\n\n\n") % Second question
 % This time we know that this function is a free falling object's movement
 % function so we write coefficients and their units
fprintf("Gravity is %f m/s^2\n\n",x(1)*-2) % Since the formula is (gt^2)/2 I multiplied with 2, 
 % also I multiplied with -1 because gravity is negative in the formula
fprintf("Starting speed of the object is %f m/s\n\n",x(2)) % Positive value indicates that object is thrown up
fprintf("Starting height of the object is %f m\n\n",x(3))
 
fprintf("The mathematical model of movement is: x(t) = ")
for ii=3:-1:2
    fprintf("%f t^%d + ",x(length(x)+1-ii),ii-1) % Height/time formula according to report   
end
fprintf("%f\n\n",x(length(x)))

interpolatedFunction=polyval(x,xi); % That part is required for plotting the calculated function

ax1=subplot(2,2,[1,2]); % I divided the plot into subplots
ax2=subplot(2,2,[3,4]);

subplot(ax1) % Comparison of measured and calculated height/time graphs
hold on
title(ax1,'Object height/time')
xlabel('Time (s)')
ylabel('Height (m)')
stem(xi,Yi,'g.','MarkerSize',5,'LineStyle','none') % Measured values are plotted here
axis([min(xi)-1 max(xi)*10/9 -5 max(Yi)*10/9]) % Limits of the graph for better visualization
plot(xi,interpolatedFunction,'r') % Calculated values are plotted here
legend('Measured Values','Interpolated Function'); % Legend of the graph


hold off
subplot(ax2) % In this graph we compare how close our calculation is to measured values
hold on
title(ax2,'Residuals')
xlabel('Time (s)')
ylabel('Deflection (m)')

bar(ax2,xi,Yi-interpolatedFunction)

deflection = sumabs(Yi-interpolatedFunction)/length(xi); % Shows the amount of deflection
normOfResiduals = norm(Yi-interpolatedFunction); % Compute the norm of the residuals 
                %(a statistic you can use to analyze how well a model fits your data)

fprintf("Average deflection from function is %f\n\n",deflection) % Comparison values are printed here
fprintf("Norm of residuals is %f\n\n\n",normOfResiduals)

hold off
enteredValue= 'Press 1 to plot the measured falling vs function falling ' ; % Lets user observe plotted graphs
% until user enters 1 and pressses enter to see the visualization of falling, velocity/time
% height/time graphs
keepGoing = input(enteredValue);


if keepGoing == 1
    close all   % Closes the existing graphs to open place for the new graphs 
end


hmax = max(interpolatedFunction); % Max height is calculated here. 
if hmax < max(interpolatedFunction) % this value is used in the visualization of graphs later on
    hmax = max(interpolatedFunction);
end
newtonFalling(hmax,xi,Yi,x)

function newtonFalling(h,xi,yi,x)

        t = 0; % Instant time of the falling simulation
        a = 1; % Integer time for measured ball height
        b = 0; % A value used to calculate integer time for measured ball height
        
        height=h; % Maximum height value for limiting the graphs
        startingHeight = x(3); % Gravity and starting values for velocity and height are defined
        startingVelocity = x(2);
        accelerationOfObject = x(1)*2;
        fallTime = max(roots(x)); % 1 negative and 1 positive roots. 
        % Since time cannot be negative % fallTime is the greater one
        
        % In order to visualise the program better I calculated max and min
        % velocity values here
        % x(t) = -4.935419 t^2 + 50.284380 t^1 + 200.144250 derivating over t is;
        % v(t) = -9.870838*t + 50.284380
        velocityLimitSolver = 0:0.01:fallTime;  
        % This part calculates the min and max speeds of the object for 
        % visualising the velocity/time graph
        velocityFunction = accelerationOfObject*velocityLimitSolver + startingVelocity;
        velocityDownLimit = min(velocityFunction)*10/9;
        velocityUpLimit = max(velocityFunction)*10/9;
        
        % 3 subplots are created for velocity/time, height/time graphs and
        % visual simulation of the falling of the objects on a graph
        ax1=subplot(2,2,1);
        ax2=subplot(2,2,2);
        ax3=subplot(2,2,[3,4]);
        while (true)
 
            hold on
            subplot(ax3)
            text(2.2,height*1/4,'Calculated','Color','green','FontSize',12)
            text(6.2,height*1/4,'Measured','Color','red','FontSize',12)
            plot(7,yi(a),'r.','MarkerSize',30) % This part plots the free falling of measured object
            
            % This part plots the free falling of calculated object
            plot(3,(accelerationOfObject*t^2)/2 + startingVelocity*t + startingHeight ,'g.','MarkerSize',30)
            
            title('Ball falling simulation');
            ylabel('Height (m)');
            xlabel('Ground');
            ylim([0 height*10/9]);% In order to visualize the falling upper 
                                  % limit of graph needs to stays still
            xlim([-0 10]);
            hold off
            
            heightOfRed(b+1) = yi(a);
            if (mod(b,10)==0) && a < length(yi) % b's increment ratio is 10 times of the increment
                % of t so as t and b increase 10 times 1 second is passed
                % This part is used to calculate location of measured ball
                a = a+1;
            end
            
            hold on
            subplot(ax1) % Height/time graph for calculated object
            title('Calculated height of the ball');
            ylabel('Height (m)');
            xlabel('Time(s)');
            xlim([-.3 fallTime*10/9]); % Limitting for visualization
            ylim([0 height*10/9]);
            c = 0:.1:t; % a is for making graphs of distance and velocity over 
                        % time. Since t is a constant number
                        % plot function wouldn't work beause of unmatched sizes
                        % so we created a matrice
            plot(c,(accelerationOfObject*c.^2)/2 + startingVelocity*c + startingHeight,'g')
            hold off                    
            
            hold on
            subplot(ax2) % Velocity/time graph for calculated object
            plot(c,startingVelocity+accelerationOfObject*c,'g'); % Derivating time in (g*t^2)/2 to find change 
                                                                 % in velocity -> g*t

            title('Calculated velocity of the ball');
            ylabel('Velocity (m/s)');
            xlabel('Time(s)');
            xlim([-.3 fallTime*10/9]); % Limitting for visualization
            ylim([velocityDownLimit velocityUpLimit]);

            hold off
                                 
                                 
            t = t + 0.1; % t and b increase
            b = b + 1; % t is used for graphs of calculated free falling
            % b is used for plotting measured value's free falling
            
            pause(0.04) % Pause for our human eyes to see the simulation
                        % and adding aspect of reality
            cla(ax3) % Clearing the objects on ax3 so they don't overlap
            if fallTime <= t % When falltime is reached loop breaks
                break
            end
        
        end
    hold on
    subplot(ax3)
    % Fallen objects are printed on the ground because of cla function they
    % need to be printed on the ground after while loop ends
    plot(ax3,3,0,'g.','MarkerSize',30);
    plot(ax3,7,0,'r.','MarkerSize',30)
    title('Ball falling simulation');
    ylabel('Height (m)');
    xlabel('Ground');
    ylim([0 height*10/9]);
    xlim([-0 10]);
    text(2.2,height*1/4,'Calculated','Color','green','FontSize',12)
    text(6.2,height*1/4,'Measured','Color','red','FontSize',12)
    
    hold off
    
    % Falling times for measured and calculated objects are printed here
    fprintf("Calculated object hit the ground after %.2f seconds\n\n",fallTime)
    fprintf("Measured object hit the ground after %.2f-%.2f seconds\n\n ",max(xi),max(xi)+1)
    

end

