function ydot=Toy06_rev_enzyme__Chain_4(t,y,params)

% parameters
% #1: reaction 1, parameter k
% #2: reaction 2, parameter Vmaxf
% #3: reaction 2, parameter KmS
% #4: reaction 2, parameter Vmaxb
% #5: reaction 2, parameter KmP
% #6: reaction 3, parameter Vmaxf
% #7: reaction 3, parameter KmS
% #8: reaction 3, parameter Vmaxb
% #9: reaction 3, parameter KmP
% #10: reaction 4, parameter Vmaxf
% #11: reaction 4, parameter KmS
% #12: reaction 4, parameter Vmaxb
% #13: reaction 4, parameter KmP
% #14: reaction 5, parameter k

% equations
ydot=zeros(5,1);
ydot(1)=+1*(params(1))-1*(((params(2)/params(3)*y(1)) - (params(4)/params(5)) * y(2)) / (1+y(1)/params(3)+y(2)/params(5)));
ydot(2)=+1*(((params(2)/params(3)*y(1)) - (params(4)/params(5)) * y(2)) / (1+y(1)/params(3)+y(2)/params(5)))-1*(((params(6)/params(7)*y(2)) - (params(8)/params(9)) * y(3)) / (1+y(2)/params(7)+y(3)/params(9)));
ydot(3)=+1*(((params(6)/params(7)*y(2)) - (params(8)/params(9)) * y(3)) / (1+y(2)/params(7)+y(3)/params(9)))-1*(((params(10)/params(11)*y(3)) - (params(12)/params(13)) * y(4)) / (1+y(3)/params(11)+y(4)/params(13)));
ydot(4)=+1*(((params(10)/params(11)*y(3)) - (params(12)/params(13)) * y(4)) / (1+y(3)/params(11)+y(4)/params(13)))-1*(params(14)*y(4));
ydot(5)=+1*(params(14)*y(4));
