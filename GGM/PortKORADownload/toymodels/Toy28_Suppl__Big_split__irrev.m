function ydot=Toy28_Suppl__Big_split__irrev(t,y,params)

% parameters
% #1: reaction 1, parameter k
% #2: reaction 2, parameter k
% #3: reaction 3, parameter k
% #4: reaction 4, parameter k
% #5: reaction 5, parameter k
% #6: reaction 6, parameter k
% #7: reaction 7, parameter k
% #8: reaction 8, parameter k

% equations
ydot=zeros(6,1);
ydot(1)=+1*(params(1))-1*(params(2)*y(1));
ydot(2)=+1*(params(2)*y(1))-1*(params(3)*y(2))-1*(params(6)*y(2));
ydot(3)=+1*(params(3)*y(2))-1*(params(4)*y(3));
ydot(4)=+1*(params(4)*y(3))-1*(params(5)*y(4));
ydot(5)=+1*(params(6)*y(2))-1*(params(7)*y(5));
ydot(6)=+1*(params(7)*y(5))-1*(params(8)*y(6));
