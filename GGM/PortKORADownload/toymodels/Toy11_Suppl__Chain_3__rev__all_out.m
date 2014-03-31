function ydot=Toy11_Suppl__Chain_3__rev__all_out(t,y,params)

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
ydot=zeros(3,1);
ydot(1)=+1*(params(1))-1*(params(2)*y(1))+1*(params(3)*y(2))-1*(params(6)*y(1));
ydot(2)=+1*(params(2)*y(1))-1*(params(3)*y(2))-1*(params(4)*y(2))+1*(params(5)*y(3))-1*(params(7)*y(2));
ydot(3)=+1*(params(4)*y(2))-1*(params(5)*y(3))-1*(params(8)*y(3));
