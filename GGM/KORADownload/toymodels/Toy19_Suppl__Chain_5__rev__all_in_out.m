function ydot=Toy19_Suppl__Chain_5__rev__all_in_out(t,y,params)

% parameters
% #1: reaction 1, parameter k
% #2: reaction 2, parameter k
% #3: reaction 3, parameter k
% #4: reaction 4, parameter k
% #5: reaction 5, parameter k
% #6: reaction 6, parameter k
% #7: reaction 7, parameter k
% #8: reaction 8, parameter k
% #9: reaction 9, parameter k
% #10: reaction 10, parameter k
% #11: reaction 11, parameter k
% #12: reaction 12, parameter k
% #13: reaction 13, parameter k
% #14: reaction 14, parameter k
% #15: reaction 15, parameter k
% #16: reaction 16, parameter k
% #17: reaction 17, parameter k
% #18: reaction 18, parameter k

% equations
ydot=zeros(5,1);
ydot(1)=+1*(params(1))-1*(params(2)*y(1))+1*(params(3)*y(2))-1*(params(15)*y(1));
ydot(2)=+1*(params(2)*y(1))-1*(params(3)*y(2))-1*(params(4)*y(2))+1*(params(5)*y(3))+1*(params(11))-1*(params(16)*y(2));
ydot(3)=+1*(params(4)*y(2))-1*(params(5)*y(3))-1*(params(6)*y(3))+1*(params(7)*y(4))+1*(params(12))-1*(params(17)*y(3));
ydot(4)=+1*(params(6)*y(3))-1*(params(7)*y(4))-1*(params(8)*y(4))+1*(params(9)*y(5))+1*(params(13))-1*(params(18)*y(4));
ydot(5)=+1*(params(8)*y(4))-1*(params(9)*y(5))-1*(params(10)*y(5))+1*(params(14));
