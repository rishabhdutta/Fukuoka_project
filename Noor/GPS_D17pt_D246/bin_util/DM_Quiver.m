function H=DM_Quiver(XY,d,dcov,scale,origin,color)
% XY is 2*n where n is number of stations , d is 3n*1 that is
% 3n datapoints, dcov is in 3n*3n matrix format
%Adds velocity vector and error ellipses to the map

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   May 23 2005   Sjonni                Changes
%   Sep 9, 2000   Peter                 Major re-write
%   April, 2000   Peter Cervelli        Original Code
%
%-------------------------------------------------------------

ns=size(XY,2);

%Arrow shafts

   switch length(d)==ns

      case 0 %Horizontal

         U=d(1:3:end)'*scale*1e2;
	 U2=d(1:3:end)'*scale*1e2/1.15;
         V=d(2:3:end)'*scale*1e2;
	 V2=d(2:3:end)'*scale*1e2/1.15;

      case 1 %Vertical

         V=d'*scale*1e2;
         V2=d'*scale*1e2/1.15;
         U=zeros(size(V)); U2=U;
   end

   X=[XY(1,:);XY(1,:)+U;repmat(NaN,1,ns)];
   Y=[XY(2,:);XY(2,:)+V;repmat(NaN,1,ns)];

%Arrow heads

   alpha=0.2;
   beta=0.33;
   Up=U;
   Vp=V;
   L=sqrt(sum([X(1,:)-X(2,:);Y(1,:)-Y(2,:)].^2));
   I=find(L>3);            % Avoid too large arrow heads
   Up(I)=Up(I)./(L(I)/3);
   Vp(I)=Vp(I)./(L(I)/3);
  
   I2=find(L<0.8);L(I2);    % Avoid too small arrow heads
   Up(I2)=Up(I2)./(L(I2)/0.8); U2(I2)=U(I2)-(U2(I2)./L(I2))*0.1;
   Vp(I2)=Vp(I2)./(L(I2)/0.8); V2(I2)=V(I2)-(V2(I2)./L(I2))*0.1;
   % Arrow and arrow head:
   X=[X;X(2,:)-alpha*(Up+beta*(Vp+eps));X(2,:);X(2,:)-alpha*(Up-beta*(Vp+eps));XY(1,:)+U2;X(2,:)-alpha*(Up+beta*(Vp+eps));repmat(NaN,1,ns)];
   Y=[Y;Y(2,:)-alpha*(Vp-beta*(Up+eps));Y(2,:);Y(2,:)-alpha*(Vp+beta*(Up+eps));XY(2,:)+V2;Y(2,:)-alpha*(Vp-beta*(Up+eps));repmat(NaN,1,ns)];
   
   
%Error ellipses

   if ~isempty(dcov)
      switch length(d)==ns

         case 0 %Horizontal

            r=linspace(0,2*pi,100);
            for i=ns:-1:1
               I=(i-1)*3+1;
               DCOV=scale^2*1e4*full(dcov(I:I+1,I:I+1));
               [v,w]=eig(DCOV);
               az=pi/2-atan2(v(1,1),v(2,1));
               w=sqrt(diag(w)*5.9915);  % Comes from 2D statistics
               elpts=[cos(az) -sin(az);sin(az) cos(az)]*[w(1)*cos(r);w(2)*sin(r)];
               X(10:110,i)=[elpts(1,:)'+X(2,i);NaN];
               Y(10:110,i)=[elpts(2,:)'+Y(2,i);NaN];
            end
            
         case 1 %Vertical, make error bars using another color
               DVAR = diag(dcov);           % get variances
	       VSTD = sqrt(DVAR(3:3:end));  % get vertical st.dev.
	       Vstd = VSTD'*scale*1e2 * 2;  % scale std values,
                                            % and take 95%
	       AX = X(2,:); AY=Y(2,:);      % get arrow ends
	       hb = 0.5*scale;
	       XV = [AX;AX;AX-hb;AX+hb;AX;AX;AX-hb;AX+hb;NaN(1,ns)];
	       YV = [AY;AY+Vstd;AY+Vstd;AY+Vstd;AY+Vstd;AY-Vstd;AY-Vstd;AY-Vstd;repmat(NaN,1,ns)];
      end

   end

%Change to LLH if necessary

 %  if nargin > 4
 %     LLH=DM_local2llh([X(:)';Y(:)'],origin);
 %     X=LLH(1,:);
 %     Y=LLH(2,:);
 %  end
   
   if nargin < 5
     color = 'k';
   end
   
   color2='c';
   
%Add vector field to map


 
   
   if ~isempty(dcov)
      switch length(d)==ns

         case 0 %Horizontal
   
         case 1 %Vertical, make error bars using another color
         
      H=line('XData',        XV(:), ...
          'YData',        YV(:), ...
          'ZData',        zeros(prod(size(XV)),1), ...
          'Tag',          'StdVector', ...
          'HitTest',      'off', ...
          'Color',        color2);
      end

   end




   H=line('XData',        X(:), ...
          'YData',        Y(:), ...
          'ZData',        zeros(prod(size(X)),1), ...
          'Tag',          'DataVector', ...
          'HitTest',      'off', ...
          'Color',        color);
   fill(X(4:8,:),Y(4:8,:),color,'EdgeColor','none');
   
   
   
