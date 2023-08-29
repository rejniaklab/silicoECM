function TACS_patterns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a companion code for paper "Dynamics of fibril collagen %%
%% remodeling by tumor cells: a model of tumor-associated collagen %%
%% signatures" by S. Poonja, A. Forero Pinto, M. Lloyd, M. Damaghi %%
%% and K.A. Rejniak                                                %% 
%% This code simulates development of tumor-associated collagen    %%
%% signatures (TACS) during the growth of a tumor cell cluster.    %%
%% It generates different ECM and tumor morphologies as shown in   %%
%% Figure 9 of the paper. The following parameters need to be      %%
%% specified:                                                      %%
%%  alpha : from [0,1] 1=persistent migration; 0=contact guidance  %% 
%%  beta  : from [0,1] 0=no ECM remodeling; 1=aligned with cell    %% 
%%  divThr: ECM stiffness level for growth arrest;                 %%  
%%  motThr: ECM stiffness for starting cell migration;             %%
%%                                                                 %%
%% August 27th, 2023                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -----------------------------%
% model parameters to specify:
alpha=1;     % in [0,1]; persistence migration coefficient
beta =0.5;   % in [0,1]; ECM compliance to remodeling coefficient 
divThr=50;   % in [10,70]; ECM stiffness for growth-arrest
motThr=25;   % in [10,150]; ECM stiffness for cell migration



% -----------------------------%
% no changes in the code below

% computational domain
xmin=-250; xmax=-xmin; ymin=xmin; ymax=xmax;   % domain boundries
dx=4; dy=4;                                    % grid width 
[xx,yy]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);   % domain for drawing
Nx=length(xmin:dx:xmax);                       % grid size
Ny=length(ymin:dy:ymax);      

% time parameters 
time=6;                    % time of simulations [days]   
dt  =2.5;                  % time step [sec]
Niter= time*24*60*60/dt;   % total number of iterations [days]

% tumor cell parameters and variables
cell=[0,0];                % cell locations [x,y]
Ncell=1;                   % number of cells
cellAge=[10*60*60];        % cell current age [sec]
matAgeInit=20*60*60;       % cell default maturation age [sec]
matAgeSTD=0.15;            % std of cell maturation age 15% 
matAge=[matAgeInit];       % individual cell maturation age [sec]
rad=8;                     % cell radius [microns]
cellmotile=zeros(Ncell,1); % 0=cell is not motile, 1=cell is motile

% tumor cell neighborhood parameters
neigh=zeros(Ncell,1);      % current number of cell neighbors
Neighrad = 4*rad;          % overcrowding neighborhood radius [microns]
neighMax=12;               % overcrowding max number of neighbors
neighEdgeMax=5;            % max number of neighbors on the edge for adhesion

decm=4;        % radius for remodeling ECM when a growing cell pushes
dpull=1;       % radius for remodeling ECM when a migrating cell pulls
dmot=2;        % radius for sensing ECM when cell is migrating
dsens=2;       % radius for sensing ECM for cell division

% cell-cell and cell-ECM interaction forces
stiffR=50;             % repulsive force stiffness [ug/(um.s^2)]
stiffA= 5;             % adhesive force stiffness [ug/(um.s^2)]
adhRmin=2.2;           % min radius for edge cell adhesion [um]
adhRmax=2*adhRmin;     % max radius for edge cell adhesion [um]
speedM= 0.08;          % motility force magnitude
stiff_inc=0.0001/dt;   % ECM stiffness increase when push/pull
nu=250;                % medium viscosity [ug/(um.s)]


% ECM initiation -- random direction, uniform stiffness
ecmStiff=ones(Nx,Ny);            % ecm stiffness [default=1]
ecmy=2*(0.5-rand(Nx,Ny));        % ecm direction [random]
ecmx=2*(0.5-rand(Nx,Ny));



%------- MAIN LOOP ---------%
for iter=0:Niter

  % report simulation progress  
  if mod(iter,6*60*60/dt)==0
    disp(['simulated time passed: ',num2str(iter*dt/3600),...
          ' hours out of ',num2str(Niter*dt/3600),' hours'])
  end


  % cell aging  
  cellAge=cellAge+dt;  
  
  
  % number of neighbors for proliferation
  neigh=zeros(Ncell,1);    
  for jj=1:Ncell
    for kk=jj+1:Ncell
      dxx=cell(jj,1)-cell(kk,1); dyy=cell(jj,2)-cell(kk,2);
      dist=sqrt(dxx^2+dyy^2);
      if (dist>0)&&(dist<=Neighrad)   
        neigh(jj,1)=neigh(jj,1)+1; neigh(kk,1)=neigh(kk,1)+1;
      end
    end
  end
  

  % ECM stiffness sensed by the cell
  cellECM=zeros(Ncell,1);
  for jj=1:Ncell
    ix=floor((cell(jj,1)-xmin)/dx);  iy=floor((cell(jj,2)-ymin)/dy);
    for ik=-dsens:dsens
      for im=-dsens:dsens
        xdist=(cell(jj,1)-(xmin+(ix+ik)*dx)); 
        ydist=(cell(jj,2)-(ymin+(iy+im)*dy));  
        dist=sqrt(xdist^2+ydist^2);       
        if (dist>rad)&&(dist<=rad+2*dx)
          if (1+iy+im>=1)&&(1+iy+im<=Nx)&&(1+ix+ik>=1)&&(1+ix+ik<=Ny) 
            cellECM(jj,1)=cellECM(jj,1)+ecmStiff(1+iy+im,1+ix+ik);
          end %if
        end %if
      end %for im
    end %for ik
  end %for jj

 
  % condition for cell becoming motile
  ind=find((cellECM>motThr)&(neigh<=neighMax)); 
  if length(ind)>0 cellmotile(ind)=1; end


  % condition for cell division
  for ii=1:Ncell         
    if (cellAge(ii)>=matAge(ii))&&(neigh(ii,1)<=neighMax)&&...
       (cellECM(ii)<=divThr)&&(cellmotile(ii)==0)     
                   
       Ncell=Ncell+1;   % new daughter cell
       th=rand*2*pi;    % daughter cells placement  
       cell(Ncell,1:2)=cell(ii,1:2)+0.5*rad*[cos(th),sin(th)];
       cell(ii,1:2)   =cell(ii,1:2)-0.5*rad*[cos(th),sin(th)];

       cellAge(Ncell)=0;  cellAge(ii)=0; % daughter cell ages
       matAge(Ncell)=matAge(ii)+matAgeSTD*(2*rand-1)*matAge(ii);
       matAge(ii)   =matAge(ii)+matAgeSTD*(2*rand-1)*matAge(ii); 

       cellmotile(Ncell)=0;  % non-motile by default
       cellECM(Ncell)=0;     % no ECM sensing by default
    end 
  end   
  

  % cell-cell repulsive interactions
  forceR=zeros(Ncell,2);   % repulsive forces
  neigh=zeros(Ncell,1);    % number of neighbors
  for jj=1:Ncell
    for kk=jj+1:Ncell
      dxx=cell(jj,1)-cell(kk,1); dyy=cell(jj,2)-cell(kk,2);
      dist=sqrt(dxx^2+dyy^2);
      if (dist>0)&&(dist<2*rad)          % repulsive forces
        forceR(jj,1:2)=forceR(jj,1:2)+stiffR*(2*rad-dist)*[dxx,dyy]/dist;
        forceR(kk,1:2)=forceR(kk,1:2)-stiffR*(2*rad-dist)*[dxx,dyy]/dist;
      end        
      if (dist>0)&&(dist<=Neighrad)      % neighborhood
        neigh(jj,1)=neigh(jj,1)+1;
        neigh(kk,1)=neigh(kk,1)+1;
      end
    end
  end

  
  % cell-cell adhesive interactions for non-motile cells on the edge
  forceA=zeros(Ncell,2);   % adhesive forces
  ind=find(neigh<=neighEdgeMax);  % only cells on the edge
  for jj=1:length(ind)
    if cellmotile(ind(jj))==0  
      for kk=1:Ncell
        if cellmotile(kk)==0  
          dxx=cell(ind(jj),1)-cell(kk,1); dyy=cell(ind(jj),2)-cell(kk,2);
          dist=sqrt(dxx^2+dyy^2);
          if (dist>adhRmin*rad)&&(dist<adhRmax*rad)    % adhesive forces
            forceA(ind(jj),1:2)=forceA(ind(jj),1:2)+stiffA*(adhRmin*rad-dist)*[dxx,dyy]/dist;
            forceA(kk     ,1:2)=forceA(kk,1:2)     -stiffA*(adhRmin*rad-dist)*[dxx,dyy]/dist;
          end % if
        end %if cellmotile
      end %for kk
    end %if cellmotile
  end %for jj
  
  
  % cell motility force
  fmotile=zeros(Ncell,2);   % motility forces
  ind=find(cellmotile==1);  % only for motile cells
  for jj=1:length(ind) 
    ix=floor((cell(ind(jj),1)-xmin)/dx);  
    iy=floor((cell(ind(jj),2)-ymin)/dy);
    for ik=-dmot:dmot
      for im=-dmot:dmot
        ixx=1+ix+ik; iyy=1+iy+im;  % sens ECM directions nearby  
        if (1+iy+im>=1)&&(1+iy+im<=Nx)&&(1+ix+ik>=1)&&(1+ix+ik<=Ny)        
          fmotile(ind(jj),1)=fmotile(ind(jj),1)+ecmx(iyy,ixx);           
          fmotile(ind(jj),2)=fmotile(ind(jj),2)+ecmy(iyy,ixx);      
        end %if
      end % for im
    end % for ik
    xdist=cell(ind(jj),1); ydist=cell(ind(jj),2); % persistent direction
    xydist=sqrt(xdist^2+ydist^2); 
    if xydist>0     % combination of ECM guided & persistent movement
      fmotile(ind(jj),1)=(1-alpha)*fmotile(ind(jj),1)+alpha*xdist;        
      fmotile(ind(jj),2)=(1-alpha)*fmotile(ind(jj),2)+alpha*ydist;   
    end
    dmx=fmotile(ind(jj),1); dmy=fmotile(ind(jj),2); 
    dmxy=sqrt(dmx*dmx+dmy*dmy);
    if dmxy>0       % force=magnitude * unit direction 
      fmotile(ind(jj),1:2)=speedM*fmotile(ind(jj),1:2)/dmxy;
    end
  end % for


  % ECM remodeling due to cell pushing [tangential orientation]
  in_ind=[]; % find all inner points
  for jj=1:Ncell
    if cellmotile(jj)==0  % non-motile cells
      ix=floor((cell(jj,1)-xmin)/dx);  iy=floor((cell(jj,2)-ymin)/dy);
      for ik=-decm:decm
        for im=-decm:decm
          xdist=(cell(jj,1)-(xmin+(ix+ik)*dx));
          ydist=(cell(jj,2)-(ymin+(iy+im)*dy));  
          dist=sqrt(xdist^2+ydist^2);       
          if (1+iy+im>=1)&&(1+iy+im<=Nx)&&(1+ix+ik>=1)&&(1+ix+ik<=Ny) 
            if dist<=rad  % find grid indices inside the cells
              in_ind=[in_ind; 1+ix+ik,1+iy+im]; 
            end %if dist
          end % if
        end % for im
      end % for ik
    end % if motile 
  end  % for jj

  
  % change ECM around the pushing or pulling cells
  for jj=1:Ncell
    ix=floor((cell(jj,1)-xmin)/dx);  iy=floor((cell(jj,2)-ymin)/dy);
    
    if cellmotile(jj)==1   % motile cells will pull ECM
      for ik=-dpull:dpull
        for im=-dpull:dpull
          ixx=1+ix+ik; iyy=1+iy+im;  
          if (1+iy+im>=1)&&(1+iy+im<=Nx)&&(1+ix+ik>=1)&&(1+ix+ik<=Ny)   
            % change direction of ECM for compliance coefficient beta  
            ecmx(iyy,ixx)=(1-beta)*ecmx(iyy,ixx)+beta*fmotile(jj,1);          
            ecmy(iyy,ixx)=(1-beta)*ecmy(iyy,ixx)+beta*fmotile(jj,2); 
              
            % change stiffness of 3 layers of ECM arounf the cell 
            xdist=cell(jj,1)-(xmin+ixx*dx); ydist=cell(jj,2)-(ymin+iyy*dx);
            dist=sqrt(xdist*xdist+ydist*ydist);
            if (dist<=rad)
              ecmStiff(iyy,ixx)=ecmStiff(iyy,ixx)+3*beta*stiff_inc;
            elseif (dist<=rad+dx*sqrt(2))
              ecmStiff(iyy,ixx)=ecmStiff(iyy,ixx)+2*beta*stiff_inc;
            elseif (dist<=rad+2*dx*sqrt(2))
              ecmStiff(iyy,ixx)=ecmStiff(iyy,ixx)+1*beta*stiff_inc;
            end % if 
          end % if
        end % for im
      end % for ik

    else
 
      % non-motile cells that are relocated and push on the ECM   
      mag=sqrt((forceR(jj,1)+forceA(jj,1))^2+(forceR(jj,2)+forceA(jj,2))^2);
      if mag>0 % cell pushes on ECM (force non-zero)
        out_ind1=[]; out_ind2=[]; out_ind3=[]; %outer points per cell
        ecm_tomove=0;    
        for ik=-decm:decm
          for im=-decm:decm
            xdist=(cell(jj,1)-(xmin+(ix+ik)*dx));
            ydist=(cell(jj,2)-(ymin+(iy+im)*dy));  
            dist=sqrt(xdist^2+ydist^2);       
            if (1+iy+im>=1)&&(1+iy+im<=Nx)&&(1+ix+ik>=1)&&(1+ix+ik<=Ny) 
              if dist<=rad
                ecm_tomove=ecm_tomove+ecmStiff(1+iy+im,1+ix+ik);
                ecmx(1+iy+im,1+ix+ik)=0;      %  no ecm inside the cells
                ecmy(1+iy+im,1+ix+ik)=0;      
                ecmStiff(1+iy+im,1+ix+ik)=0;
              elseif dist<=rad+dx*sqrt(2)     % 3 layers of ECM near cell
                out_ind1=[out_ind1; 1+ix+ik,1+iy+im];
              elseif dist<=rad+dx+dx*sqrt(2)
                out_ind2=[out_ind2; 1+ix+ik,1+iy+im];
              elseif dist<=rad+2*dx+dx*sqrt(2)
                out_ind3=[out_ind3; 1+ix+ik,1+iy+im];                
              end %if dist
            end % if
          end % for im
        end % for ik
 
        % determine if the grid points in out_ind are not inside the cells
        out_ind1=index_diff(in_ind  ,out_ind1);
        out_ind2=index_diff(in_ind  ,out_ind2);
        out_ind2=index_diff(out_ind1,out_ind2);
        out_ind3=index_diff(in_ind,out_ind3);
        out_ind3=index_diff(out_ind1,out_ind3);
        out_ind3=index_diff(out_ind2,out_ind3);
        Nfiber=3*size(out_ind1,1)+2*size(out_ind2,1)+size(out_ind3,1);
  
        if Nfiber>0    % if there are fibers to push from cell location
          ecm_tomove=ecm_tomove/Nfiber;
   
          for ii=1:size(out_ind1,1)
            ixx=out_ind1(ii,1); iyy=out_ind1(ii,2);   
            xdist=(cell(jj,1)-(xmin+ixx*dx));
            ydist=(cell(jj,2)-(ymin+iyy*dy));  
            % change direction of ECM for compliance coefficient beta 
            ecmx(iyy,ixx)=beta*ecmx(iyy,ixx)+(1-beta)*ydist;           
            ecmy(iyy,ixx)=beta*ecmy(iyy,ixx)+(1-beta)*(-1)*xdist;  
            % change ECM stiffness
            ecmStiff(iyy,ixx)=ecmStiff(iyy,ixx)+3*ecm_tomove+beta*stiff_inc;  
          end
          for ii=1:size(out_ind2,1)
            ixx=out_ind2(ii,1); iyy=out_ind2(ii,2);   
            xdist=(cell(jj,1)-(xmin+ixx*dx));
            ydist=(cell(jj,2)-(ymin+iyy*dy));  
            % change direction of ECM for compliance coefficient beta 
            ecmx(iyy,ixx)=beta*ecmx(iyy,ixx)+(1-beta)*ydist;           
            ecmy(iyy,ixx)=beta*ecmy(iyy,ixx)+(1-beta)*(-1)*xdist;   
            % change ECM stiffness
            ecmStiff(iyy,ixx)=ecmStiff(iyy,ixx)+2*ecm_tomove+beta*stiff_inc;  
          end
          for ii=1:size(out_ind3,1)
            ixx=out_ind3(ii,1); iyy=out_ind3(ii,2);
            xdist=(cell(jj,1)-(xmin+ixx*dx));
            ydist=(cell(jj,2)-(ymin+iyy*dy));  
             % change direction of ECM for compliance coefficient beta 
            ecmx(iyy,ixx)=beta*ecmx(iyy,ixx)+(1-beta)*ydist;           
            ecmy(iyy,ixx)=beta*ecmy(iyy,ixx)+(1-beta)*(-1)*xdist;    
            % change ECM stiffness
            ecmStiff(iyy,ixx)=ecmStiff(iyy,ixx)+ecm_tomove+beta*stiff_inc; 
          end  
        end % if Nfiber>0
      end % if mag
    end % if-else cellmotile
  end  % for jj
  
  % define ECM direction as a unit vector (magnitude is in ecmStiff)
  magnitude=sqrt(ecmx.^2+ecmy.^2);
  for ia=1:Nx
    for ib=1:Ny
      if (magnitude(ia,ib)>0)
        ecmx(ia,ib)=ecmx(ia,ib)/magnitude(ia,ib);
        ecmy(ia,ib)=ecmy(ia,ib)/magnitude(ia,ib);
      end
    end
  end
  
  
  % cell relocation due to cummulative force: repulsive+adhesive+motility
  forceTot=forceR+forceA+fmotile; 
  cell=cell+(forceTot)*dt/nu;

  
end %main loop
% ----- main loop end ----- %



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ----- draw cell cluster and ECM after simulation  

disp('  ~~~ the figure is being generated')

ff=figure('position',[500,100,650,650]);  % make figure
set(ff,'PaperPositionMode','auto');
set(ff,'InvertHardcopy','off');
clf
drawnow

% color map for drawing fiber stiffnes
mapHE=[1,0.8,1;1,0.8,1;1,0.6,1;1,0.6,1;1,0.4,1;1,0.4,1;...
       1,0,1;0.8,0,0.8;0.6,0,0.6;0.4,0,0.4];  

% fill the domain in white
fill([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],'w');
hold on  
    
  % define length of the ECM vectors for drawing with curvvec
  ecmx=ecmx.*ecmStiff; ecmy=ecmy.*ecmStiff;
  curvvec(xx,yy,ecmx,ecmy,'color',mapHE,'minmag',0.25,'maxmag',2.5)
  colormap(mapHE);    % color map
  colorbar
  caxis([1,round(max(max(ecmStiff)))+2])
    
  % define 3D cells with light on
  [sx,sy,sz]=sphere(25);
  for ii=1:Ncell  
    xp=cell(ii,1)+1.1*rad*sx; yp=cell(ii,2)+1.1*rad*sy; zp=0+1.2*rad*sz;
    surf(xp,yp,zp,'facecolor','m','edgecolor','none')
    hold on
  end
  axis([xmin,xmax,ymin,ymax,-2*rad,2*rad])
  view(3)
  light('position',[-500,500,500],'style','infinite')
  light('position',[-200,-200,200],'style','infinite')

  % flat 2D view of the tissue
  view(2)
  axis([xmin,xmax,ymin,ymax])
  axis equal
  axis([xmin,xmax,ymin,ymax])

  title([num2str(iter*dt/(24*60*60)),' days    Ncell=',num2str(Ncell),...
    ';  ECM=[',num2str(min(min(ecmStiff(3:Nx-3,3:Ny-3)))),';',...
    num2str(max(max(ecmStiff(3:Nx-3,3:Ny-3)))),']  \alpha=',...
    num2str(alpha),'   \beta=',num2str(beta),'  divThr=',...
    num2str(divThr),'  motThr=',num2str(motThr)])
  pause(0.1)
    
end % function TACS_pattern


%--------------------------------------------------
% remove from set ind_out all indices in set ind_in 
%--------------------------------------------------
function out_set=index_diff(ind_in,ind_out)
  out_set=[];
 
  for ii=1:size(ind_out,1)
    ix=ind_out(ii,1); iy=ind_out(ii,2);
    answ=0;
    for jj=1:size(ind_in,1)
      if (ix==ind_in(jj,1))&&(iy==ind_in(jj,2))
        answ=1;
      end
    end
    if answ==0
      out_set=[out_set;ix,iy];  
    end
  end
end
%--------------------------------------------------




