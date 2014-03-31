function ydot=Toy07_rev_enzyme__Chain_4__with_inhibition__Figure_1E_(t,y,params)

% parameters
% #1: reaction 1, parameter k
% #2: reaction 2, parameter Vmaxf
% #3: reaction 2, parameter KmS
% #4: reaction 2, parameter Vmaxb
% #5: reaction 2, parameter KmP
% #6: reaction 2, parameter Ki
% #7: reaction 2, parameter Kii
% #8: reaction 3, parameter Vmaxf
% #9: reaction 3, parameter KmS
% #10: reaction 3, parameter Vmaxb
% #11: reaction 3, parameter KmP
% #12: reaction 4, parameter Vmaxf
% #13: reaction 4, parameter KmS
% #14: reaction 4, parameter Vmaxb
% #15: reaction 4, parameter KmP
% #16: reaction 5, parameter k

% equations
ydot=zeros(4,1);
ydot(1)=+1*(params(1))-1*(((params(2)/params(3)*y(1)) - (params(4)/params(5)) * y(2)) / (1+y(3)/params(6)+(y(1)/params(3)+y(2)/params(5))*(1+(y(3)/params(7)))));
ydot(2)=+1*(((params(2)/params(3)*y(1)) - (params(4)/params(5)) * y(2)) / (1+y(3)/params(6)+(y(1)/params(3)+y(2)/params(5))*(1+(y(3)/params(7)))))-1*(((params(8)/params(9)*y(2)) - (params(10)/params(11)) * y(4)) / (1+y(2)/params(9)+y(4)/params(11)));
ydot(3)=+1*(((params(12)/params(13)*y(4)) - (params(14)/params(15)) * y(3)) / (1+y(4)/params(13)+y(3)/params(15)))-1*(params(16)*y(3));
ydot(4)=+1*(((params(8)/params(9)*y(2)) - (params(10)/params(11)) * y(4)) / (1+y(2)/params(9)+y(4)/params(11)))-1*(((params(12)/params(13)*y(4)) - (params(14)/params(15)) * y(3)) / (1+y(4)/params(13)+y(3)/params(15)));
