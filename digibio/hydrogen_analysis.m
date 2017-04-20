% Let D, A, H, B be the donor, the acceptor, the hydrogen, and the
% antecedent.
% 1. |D - A| < 3.5 Ang
% 2. |H - A| < 2.5 Ang
% 3. Angle(DHA) > 90 Deg
% 4. Angle(DAB) > 90 Deg
% 5. Angle(HAB) > 90 Deg
%Criteria = [3.5, 2.5, 90, 90, 90];

function [ result ] = hydrogen_analysis( D, A, H, B )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



checkOne = norm(D-A) < 3.5;
%display(norm(D-A));

checkTwo = norm(H-A) < 2.5;
%display(norm(H-A));

checkThree = radtodeg(atan2(norm(cross(D-H,A-H)),dot(D-H,A-H))) > 90;
checkFour = radtodeg(atan2(norm(cross(D-A,B-A)),dot(D-A,B-A))) > 90;
checkFive = radtodeg(atan2(norm(cross(H-A,B-A)),dot(H-A,B-A))) > 90;

%display([checkOne checkTwo checkThree checkFour checkFive]);

result = checkOne & checkTwo & checkThree & checkFour & checkFive;

end


