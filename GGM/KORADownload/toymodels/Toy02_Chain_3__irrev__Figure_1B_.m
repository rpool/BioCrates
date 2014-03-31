function ydot=Toy02_Chain_3__irrev__Figure_1B_(t,y,params)

% parameters
% #1: reaction 1, parameter k
% #2: reaction 2, parameter k
% #3: reaction 3, parameter k
% #4: reaction 4, parameter k

% equations
ydot=zeros(3,1);
ydot(1)=+1*(params(1))-1*(params(2)*y(1));
ydot(2)=+1*(params(2)*y(1))-1*(params(3)*y(2));
ydot(3)=+1*(params(3)*y(2))-1*(params(4)*y(3));
