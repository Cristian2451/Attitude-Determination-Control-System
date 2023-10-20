function output = cp_op(input)
% Cross Product Operator
output = [0, -input(3), input(2); input(3), 0, -input(1); -input(2), input(1), 0]; 
end
