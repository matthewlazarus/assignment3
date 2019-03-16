%% Assignment 3  
% Matthew Lazarus 100962142

%% Question 2: Extension of Assignment 2
% In this question, a grid is set up, with the left sided boundary set to 1V.
% Two bottle-necks of conductivities different than the rest of the grid
% are introduced. Using the maxtrix form of the problem (GV=F), the electrostatic potential
% within the region is found and plotted. Additionally, the conducitivity,
% electric field and currents are plotted as well.

% Clear all previous variables, figures, etc, to ensure that the workspace
% is clean. 
clear all
clearvars
clearvars -GLOBAL
close all

for repeat = 1:1
    %Set the length and width of the grid. 
    L=160;
    W=120;
    %b) - Vary the size of the grid
    if(repeat>1 && repeat <7)
        L = (repeat-1)*40;
        W = (repeat-1)*30;
    end
    Lb = L/5;
    Wb = W/5;
    
    %Initialize the G,B and conductivity matrices. 
    G = sparse(L*W,L*W);
    B=zeros(L*W,1);
    B(1:W,1)=1;
    condMap = zeros(W,L);

    %Populate the conductivity matrix. 
    for lCount = 1:L
        for wCount = 1:W
            if(lCount < L/2 + Lb/2 && lCount > L/2 - Lb/2 &&...
                    (wCount > W-Wb || wCount < Wb))
                condMap( wCount,lCount) = 10^-2;               
            else
                condMap(wCount,lCount) = 1;
            end
        end
    end

    % Set the diagonal of the G matrix to 1. This value will be overwritten
    % later if it is not a boundary condition. 
    for count = 1:L*W
        G(count,count)=1;
    end

    %%
    % Loop through rows and columns, if not a boundary case, set the 
    % gradient based on the sum of adjacent conductivies. 

    for col = 1:L
        if(col~=1 &&col~=L)
            for row = 1:W  
                n = row + (col -1)*W;
                if(count~=1 && row ~=1 && row~=W)                
                    rxBefore = (condMap(row,col) + condMap(row,col-1))/2.0;
                    rxAfter = (condMap(row,col) + condMap(row,col+1))/2.0;
                    ryBefore = (condMap(row,col) + condMap(row-1,col))/2.0;
                    ryAfter = (condMap(row,col) + condMap(row+1,col))/2.0;

                    nyBefore = n-1;
                    nyAfter = n+1;
                    nxBefore = row+(col-2)*W;
                    nxAfter = row+col*W;
                    G(n,n) = -(rxBefore+rxAfter+ryBefore+ryAfter); 
                    G(n, nyBefore) =ryBefore;
                    G(n, nyAfter)=ryAfter;
                    G(n, nxBefore)=rxBefore;
                    G(n, nxAfter) =rxAfter; 

                elseif(row==1)
                    %Special Case: Bottom of Grid
                    rxBefore = (condMap(row,col) + condMap(row,col-1))/2.0;
                    rxAfter = (condMap(row,col) + condMap(row,col+1))/2.0;
                    ryAfter = (condMap(row,col) + condMap(row+1,col))/2.0;

                    nyAfter = n+1;
                    nxBefore = row+(col-2)*W;
                    nxAfter = row+col*W;
                    G(n,n) = -(rxBefore+rxAfter+ryAfter); 
                    G(n, nyAfter)=ryAfter;
                    G(n, nxBefore)=rxBefore;
                    G(n, nxAfter) =rxAfter; 
                    
                elseif(row==W)
                    %Special Case: Top of Grid
                    rxBefore = (condMap(row,col) + condMap(row,col-1))/2.0;
                    rxAfter = (condMap(row,col) + condMap(row,col+1))/2.0;
                    ryBefore = (condMap(row,col) + condMap(row-1,col))/2.0;

                    nyBefore = n-1;
                    nxBefore = row+(col-2)*W;
                    nxAfter = row+col*W;
                    G(n,n) = -(rxBefore+rxAfter+ryBefore); 
                    G(n, nyBefore) =ryBefore;
                    G(n, nxBefore)=rxBefore;
                    G(n, nxAfter) =rxAfter; 
                end        
            end
        end
    end 
    V=G\B;

    %Map the voltage to original grid.
    voltMap = zeros(W,L);
    for cols = 1:L
        for rows = 1:W
            n= rows+(cols-1)*W;
            voltMap(rows,cols)=V(n);        
        end
    end
    
    %%
    % Find the Electric field knowing $E=- \nabla V$
    
    [Ex, Ey]=gradient(voltMap);
    Ex=-Ex;
    Ey=-Ey;
    
    %%
    % Find the current density knowing $J=\sigma E$
    
    Jx = condMap.*Ex;
    Jy = condMap.*Ey;
    
    %Sum the currents at both contacts (edges), take the average to find
    %the total. 
    current1 = sum(Jx(:,1));
    current2 = sum(Jx(:,L));
    totalCurrent = (current1+current2)/2;
    
    %% Part A
    % As calculated above, the total current through the contacts for L=160, W=120, Lb
    % = 32 and Wb = 24 is 0.6235. Additionally, there is no observed
    % difference between the currents at each contact. 
    if(repeat==1)    
               
        %%
        % The plot of the electric potential can be seen in the figure below.
        % The left contact is set to V=1, and there is an almost linear
        % decrease up to the right contact, which is at V=0. This linearity is
        % slightly disturbed due to the two bottle-neck regions. 
        
        figure
        surf(1:L,1:W,voltMap)
        xlabel('x')
        ylabel('y')
        zlabel('Voltage')
        title('Electric Potential of Grid')
        colorbar;
        %view(2)

        %%
        % The electric field can be seen in the figure below. The electric
        % field is strongest in the bottle-neck regions.
        
        figure;
        x = 1:L;
        y=1:W;
        quiver(x,y,Ex,Ey);
        xlabel('x')
        ylabel('y')
        title('Electric Field of Rectangular Region')       
        
    end
end


