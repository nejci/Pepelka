% Unit test

% Toy example from paper
E = [1 1 2 2 3 3; 3 3 2 2 1 1; 1 1 1 2 2 2; 1 1 nan nan nan 2]';

[labelsCons, Kcons, time] = diclens(E);