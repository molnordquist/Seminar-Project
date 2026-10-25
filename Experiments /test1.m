%%Creating function files 
for k=1:n/p
  for i=1:p
    t(i+1)=t(i)+h;
    y(i+1,:)=eulerstep(t(i),y(i,:),h);
  end
  y(1,:)=y(p+1,:);t(1)=t(p+1);
  clearpoints(head);addpoints(head,y(1,1),y(1,3))
  addpoints(tail,y(1,1),y(1,3))
  drawnow;
end
%% Euler method 
function y=eulerstep(t,x,h)
%one step of the Euler method
y=x+h*ydot(t,x);
function y = trapstep(t,x,h)
%one step of the Trapezoid Method
z1=ydot(t,x);
g=x+h*z1;
z2=ydot(t+h,g);
y=x+h*(z1+z2)/2;