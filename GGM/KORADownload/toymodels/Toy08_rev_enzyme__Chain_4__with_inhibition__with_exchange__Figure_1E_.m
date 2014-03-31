function ydot=Toy08_rev_enzyme__Chain_4__with_inhibition__with_exchange__Figure_1E_(t,y,params)

% parameters
% #1: reaction 1, parameter k
% #2: reaction 2, parameter k
% #3: reaction 3, parameter k
% #4: reaction 4, parameter k
% #5: reaction 5, parameter Vmaxf
% #6: reaction 5, parameter KmS
% #7: reaction 5, parameter Vmaxb
% #8: reaction 5, parameter KmP
% #9: reaction 6, parameter k
% #10: reaction 7, parameter k
% #11: reaction 8, parameter Vmaxf
% #12: reaction 8, parameter KmS
% #13: reaction 8, parameter Vmaxb
% #14: reaction 8, parameter KmP
% #15: reaction 9, parameter k
% #16: reaction 10, parameter k
% #17: reaction 11, parameter Vmaxf
% #18: reaction 11, parameter KmS
% #19: reaction 11, parameter Vmaxb
% #20: reaction 11, parameter KmP
% #21: reaction 11, parameter Ki
% #22: reaction 11, parameter Kii

% equations
ydot=zeros(4,1);
ydot(1)=+1*(params(1))-1*(params(2)*y(1))-1*(((params(17)/params(18)*y(1)) - (params(19)/params(20)) * y(2)) / (1+y(4)/params(21)+(y(1)/params(18)+y(2)/params(20))*(1+(y(4)/params(22)))));
ydot(2)=+1*(params(3))-1*(params(4)*y(2))-1*(((params(5)/params(6)*y(2)) - (params(7)/params(8)) * y(3)) / (1+y(2)/params(6)+y(3)/params(8)))+1*(((params(17)/params(18)*y(1)) - (params(19)/params(20)) * y(2)) / (1+y(4)/params(21)+(y(1)/params(18)+y(2)/params(20))*(1+(y(4)/params(22)))));
ydot(3)=+1*(((params(5)/params(6)*y(2)) - (params(7)/params(8)) * y(3)) / (1+y(2)/params(6)+y(3)/params(8)))+1*(params(9))-1*(params(10)*y(3))-1*(((params(11)/params(12)*y(3)) - (params(13)/params(14)) * y(4)) / (1+y(3)/params(12)+y(4)/params(14)));
ydot(4)=+1*(((params(11)/params(12)*y(3)) - (params(13)/params(14)) * y(4)) / (1+y(3)/params(12)+y(4)/params(14)))-1*(params(15)*y(4))+1*(params(16));
